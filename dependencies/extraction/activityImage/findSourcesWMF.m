function [seeds, info] = findSourcesWMF(actT, params, excludeMask)
%FINDSOURCESWMF Source seeds from the trials' aligned activity images.
%
%   [SEEDS, INFO] = FINDSOURCESWMF(ACTT, PARAMS) finds sources in ACTT
%   (rows x cols x trials), the trials' activity images in the common frame
%   (NaN where not imaged). SEEDS is N x 2 [row col], 1-indexed and subpixel,
%   strongest first. Uses params.wmfPeakth, params.wmfDip,
%   params.wmfDeblendth and params.nanThresh.
%
%   FINDSOURCESWMF(ACTT, PARAMS, EXCLUDEMASK) seeds nothing where the
%   rows x cols logical EXCLUDEMASK is true (e.g. somata), as getActImPeaks
%   does with its exclusion mask.
%
%   1. Each trial's image minus its local 5x5 median, as a robust z-score
%      (median and MAD over the image), as each trial's own peak finding does.
%   2. The combined image Z is the mean of the trials' z images, again as a
%      robust z-score; pixels imaged in fewer than 1 - nanThresh of the trials
%      are left out.
%   3. F is the whitened matched filter of Z: the detection statistic for a
%      source of the empirical shape K in noise of the empirical
%      autocorrelation R. K is the average of patches around strong, isolated
%      maxima of Z, symmetric only under a half turn, so it keeps different
%      extents along rows and columns. With too few such maxima, K is a
%      Gaussian whose widths along rows and columns are those of the strong
%      maxima of Z (from the curvature of log Z), or 1 px with too few of
%      those. R is the autocorrelation of Z away from activity, which carries
%      the correlation of the noise along scan lines. F is Z correlated with
%      R^-1 K, again as a robust z-score.
%   4. Local maxima of F above params.wmfPeakth are candidates, except those
%      in EXCLUDEMASK or active (as in 7) in fewer than MINTRIALS trials,
%      which could otherwise suppress a neighbour and then be dropped
%      themselves. From the strongest down, a candidate within MERGEPX of a stronger kept one is
%      merged into it unless F dips by at least params.wmfDip between them, and
%      by params.wmfDip/sqrt(2) in each half of the trials (alternate trials,
%      filtered alike). Two sources are kept apart only where they are resolved
%      reproducibly, so one source is not split into seeds at its edges.
%   5. A seed's position is the centroid of Z over the 3x3 pixels around its
%      maximum.
%   6. Deblending. Z is fitted as a sum of copies of K, one per seed, each
%      with its own amplitude and subpixel position (within MAXSHIFT px of the
%      seed). The strongest maximum of the matched filter of the residual
%      within MERGEPX of a seed and at least MINSEP from every seed becomes a
%      new seed if it exceeds params.wmfDeblendth (at least params.wmfPeakth),
%      and params.wmfDeblendth/sqrt(2) in the residual of each half of the
%      trials, fitted at the same positions. Refit and repeat. This resolves
%      sources too close to show a dip in F, where they reproduce in both
%      halves. No seed is added in EXCLUDEMASK. Seeds then sit at their
%      fitted positions. params.wmfDeblendth = Inf skips this step.
%   7. A seed must lie outside EXCLUDEMASK and show activity (its trial's z image smoothed with a 1 px
%      Gaussian, above ZPRES within 1 px) in at least MINTRIALS trials.
%
%   With fewer than 2*MINHALF imaged trials (e.g. one long recording), the
%   halves are too noisy to test and a source may be active in only some of
%   so few trials: two maxima need only the dip in F to stay apart, there is
%   no deblending, and a seed need only be active in one trial. SILo seeds
%   such recordings with getActImPeaks instead.
%
%   INFO: Z, F, K, R, and per seed the value of F (or of the residual's F, for
%   seeds added by deblending; value), the number of trials it is active in
%   (nTrials) and whether deblending added it (added).

MERGEPX = 6; ZPRES = 4; MINTRIALS = 2; MINHALF = 2;
[H, W, T] = size(actT);
if nargin < 3 || isempty(excludeMask)
    excludeMask = false(H, W);
end
seeds = zeros(0, 2);
info = struct('Z', [], 'F', [], 'K', [], 'R', [], 'value', zeros(0, 1), 'nTrials', zeros(0, 1), 'added', false(0, 1));
if T == 0 || all(isnan(actT(:)))
    return
end

zT = nan(H, W, T); zA = -inf(H, W, T);
g = fspecial('gaussian', 7, 1);
for t = 1:T
    a = double(actT(:, :, t));
    if all(isnan(a(:))), continue, end
    zT(:, :, t) = robustZ(a - nanmedfilt2(a, [5 5]));
    z = robustZ(filtNan(zT(:, :, t), g)); z(isnan(z)) = -inf;
    zA(:, :, t) = ordfilt2(z, 9, ones(3)); % activity within 1 px
end
nActive = @(ix) sum(zA(ix + (0:T-1)*H*W) > ZPRES, 2);
imaged = find(squeeze(any(~isnan(zT), [1 2])))';
useHalves = numel(imaged) >= 2*MINHALF;
minTrials = MINTRIALS;
if ~useHalves
    minTrials = 1;
end
Z = combineTrials(zT, params.nanThresh);
[F, K, R, flt] = wmfImage(Z);
if useHalves
    Zh = {combineTrials(zT(:, :, imaged(1:2:end)), params.nanThresh), ...
          combineTrials(zT(:, :, imaged(2:2:end)), params.nanThresh)};
    Fh = {wmfImage(Zh{1}), wmfImage(Zh{2})};
end

E = F; E(isnan(E)) = -inf;
[r, c] = find(E == ordfilt2(E, 9, ones(3)) & E > params.wmfPeakth);
ix = sub2ind([H W], r, c);
ok = ~excludeMask(ix) & nActive(ix) >= minTrials;
r = r(ok); c = c(ok); v = E(ix(ok));
[v, o] = sort(v, 'descend'); r = r(o); c = c(o);
keep = true(size(r));
for i = 2:numel(r)
    for j = find(keep(1:i-1))'
        d = hypot(r(i) - r(j), c(i) - c(j));
        if d > MERGEPX, continue, end
        tt = linspace(0, 1, max(3, ceil(4*d)));
        xs = c(i) + tt*(c(j) - c(i)); ys = r(i) + tt*(r(j) - r(i));
        resolved = v(i) - min(interp2(F, xs, ys)) >= params.wmfDip;
        for h = 1:2*useHalves
            ln = interp2(Fh{h}, xs, ys);
            resolved = resolved && min(ln(1), ln(end)) - min(ln) >= params.wmfDip/sqrt(2);
        end
        if ~resolved
            keep(i) = false;
            break
        end
    end
end
r = r(keep); c = c(keep); v = v(keep);

seeds = zeros(numel(r), 2);
for k = 1:numel(r)
    rr = max(r(k)-1, 1):min(r(k)+1, H); cc = max(c(k)-1, 1):min(c(k)+1, W);
    w = Z(rr, cc); w(isnan(w)) = 0; w = max(w, 0);
    [CC, RR] = meshgrid(cc, rr);
    if sum(w(:)) > 0
        seeds(k, :) = [sum(RR(:).*w(:)), sum(CC(:).*w(:))] / sum(w(:));
    else
        seeds(k, :) = [r(k) c(k)];
    end
end

added = false(numel(v), 1);
if useHalves && isfinite(params.wmfDeblendth) && ~isempty(seeds)
    [seeds, v, added] = deblend(Z, Zh, seeds, v, K, flt, max(params.wmfDeblendth, params.wmfPeakth), MERGEPX, excludeMask);
end

ix = sub2ind([H W], round(seeds(:, 1)), round(seeds(:, 2)));
nTrials = nActive(ix);
ok = nTrials >= minTrials & ~excludeMask(ix);
seeds = seeds(ok, :);
info = struct('Z', Z, 'F', F, 'K', K, 'R', R, 'value', v(ok), 'nTrials', nTrials(ok), 'added', added(ok));
end


function [P, v, added] = deblend(Z, Zh, P0, v, K, flt, th, mergePx, excludeMask)
% fit copies of K at the seeds and add the significant maxima of the matched
% filter of the residual near them, while they reproduce in both halves
MAXSHIFT = 1.5; MINSEP = 1;
[H, W] = size(Z); hw = (size(K, 1) - 1)/2;
kg = griddedInterpolant({-hw:hw, -hw:hw}, K, 'cubic', 'none');
valid = ~isnan(Z);
added = false(numel(v), 1);
[P, ~, Zres] = fitSeeds(Z, P0, P0, kg, hw, valid, MAXSHIFT);
fltH = {halfNorm(flt, Zh{1}), halfNorm(flt, Zh{2})};
[XX, YY] = meshgrid(1:W, 1:H);
rejected = false(H, W);
for it = 1:max(20, numel(v))
    Fr = applyFilt(flt, Zres);
    Fr(isnan(Fr)) = -inf;
    Fr(Fr < ordfilt2(Fr, 9, ones(3))) = -inf;
    dmin = inf(H, W);
    for k = 1:size(P, 1)
        dmin = min(dmin, hypot(YY - P(k, 1), XX - P(k, 2)));
    end
    Fr(dmin > mergePx | dmin < MINSEP | rejected | excludeMask) = -inf;
    [vm, ix] = max(Fr(:));
    if vm < th
        break
    end
    [pr, pc] = ind2sub([H W], ix);
    near = max(pr-1, 1):min(pr+1, H); nearC = max(pc-1, 1):min(pc+1, W);
    ok = true;
    for h = 1:2
        ok = ok && peakNear(applyFilt(fltH{h}, residualFixed(Zh{h}, P, kg, hw)), [pr pc]) >= th/sqrt(2);
    end
    if ~ok
        rejected(near, nearC) = true;
        continue
    end
    Pn = [P; pr pc]; P0n = [P0; pr pc];
    [Pn, amp, Zresn] = fitSeeds(Z, Pn, P0n, kg, hw, valid, MAXSHIFT);
    if amp(end) <= 0 || any(hypot(Pn(1:end-1, 1) - Pn(end, 1), Pn(1:end-1, 2) - Pn(end, 2)) < MINSEP)
        rejected(near, nearC) = true;
        continue
    end
    P = Pn; P0 = P0n; Zres = Zresn; v(end + 1, 1) = vm; added(end + 1, 1) = true; %#ok<AGROW>
end
P(:, 1) = min(max(P(:, 1), 1), H); P(:, 2) = min(max(P(:, 2), 1), W);
[v, o] = sort(v, 'descend'); P = P(o, :); added = added(o);
end


function [P, amp, Zres] = fitSeeds(Z, P, P0, kg, hw, valid, maxShift)
% amplitudes and subpixel positions of copies of K at the seeds, least squares
% on Z (Gauss-Newton on K and its gradient, steps of at most 0.5 px, positions
% within maxShift of P0); residual Z minus the model
[H, W] = size(Z);
for it = 1:6
    [A, pix] = design(P, kg, hw, H, W, valid, true);
    x = ridgeLS(A, Z(pix));
    a = x(1:3:end); d = -[x(2:3:end), x(3:3:end)] ./ a;
    d(a <= 0, :) = 0; d = max(min(d, 0.5), -0.5);
    P = P0 + max(min(P + d - P0, maxShift), -maxShift);
    if max(abs(d(:))) < 0.02, break, end
end
[A, pix] = design(P, kg, hw, H, W, valid, false);
amp = ridgeLS(A, Z(pix));
Zres = Z; Zres(pix) = Z(pix) - A*amp;
end


function Zr = residualFixed(Z, P, kg, hw)
% residual of Z after fitting the amplitudes of K at fixed positions P
[H, W] = size(Z); valid = ~isnan(Z);
[A, pix] = design(P, kg, hw, H, W, valid, false);
Zr = Z; Zr(pix) = Z(pix) - A*ridgeLS(A, Z(pix));
end


function x = ridgeLS(A, z)
G = A'*A; lam = 1e-6*max(full(diag(G)));
x = (G + lam*speye(size(G, 1))) \ (A'*z);
end


function [A, pix] = design(P, kg, hw, H, W, valid, deriv)
% sparse design matrix: per seed a copy of K at its subpixel position (and,
% with deriv, the gradient of K there), over the imaged pixels within the
% copies' support
N = size(P, 1); nb = 1 + 2*deriv; e = 0.25;
[oy, ox] = ndgrid(-hw:hw); oy = oy(:); ox = ox(:);
I = cell(N, 1); J = I; V = I;
for k = 1:N
    ry = round(P(k, 1)) + oy; rx = round(P(k, 2)) + ox;
    in = ry >= 1 & ry <= H & rx >= 1 & rx <= W; ry = ry(in); rx = rx(in);
    li = sub2ind([H W], ry, rx); in = valid(li); li = li(in); ry = ry(in); rx = rx(in);
    dy = ry - P(k, 1); dx = rx - P(k, 2);
    kv = nanToZero(kg(dy, dx));
    if deriv
        ky = (nanToZero(kg(dy + e, dx)) - nanToZero(kg(dy - e, dx)))/(2*e);
        kx = (nanToZero(kg(dy, dx + e)) - nanToZero(kg(dy, dx - e)))/(2*e);
        I{k} = [li; li; li];
        J{k} = (k-1)*nb + [ones(numel(li), 1); 2*ones(numel(li), 1); 3*ones(numel(li), 1)];
        V{k} = [kv; ky; kx];
    else
        I{k} = li; J{k} = repmat(k, numel(li), 1); V{k} = kv;
    end
end
I = vertcat(I{:}); J = vertcat(J{:}); V = vertcat(V{:});
pix = unique(I); map = zeros(H*W, 1); map(pix) = 1:numel(pix);
A = sparse(map(I), J, V, numel(pix), N*nb);
end


function x = nanToZero(x)
x(isnan(x)) = 0;
end


function val = peakNear(F, p)
% maximum of F over the 3x3 pixels around p
[H, W] = size(F); r = round(p(1)); c = round(p(2));
w = F(max(r-1, 1):min(r+1, H), max(c-1, 1):min(c+1, W));
val = max(w(:));
if isempty(val) || isnan(val), val = -inf; end
end


function fl = halfNorm(flt, Zh)
% the full data's filter, as a robust z-score on a half's combined image
fl = flt; fl.med = 0; fl.sc = 1;
Fr = applyFilt(fl, Zh); x = Fr(~isnan(Fr));
fl.med = median(x); fl.sc = mad(x, 1)/0.6741891400433162;
end


function F = applyFilt(flt, Z)
% the whitened matched filter flt applied to Z (NaN where not used)
[H, W] = size(Z); m = ~isnan(Z); Z0 = Z; Z0(~m) = 0;
F = real(ifft2(fft2(Z0, flt.nf(1), flt.nf(2)) .* flt.KS));
F = (F(1:H, 1:W) - flt.med)/flt.sc; F(~m) = nan;
end


function [F, K, R, flt] = wmfImage(Z)
% whitened matched filter of the combined z image Z (NaN where not used); flt
% applies the same filter and scaling to other images
KHW = 4; AHW = 7; CLIP = 3; SFLOOR = 0.02;
[H, W] = size(Z); m = ~isnan(Z);
K = sourceKernel(Z, KHW);
% noise autocorrelation: Z away from activity (pixels above CLIP, dilated, left out)
act = imdilate(Z > CLIP, strel('disk', 3));
M = double(m & ~act); B = Z; B(M == 0) = 0;
nf = 2.^nextpow2([H W] + AHW + 1);
acf = real(ifft2(abs(fft2(B, nf(1), nf(2))).^2)) ./ max(real(ifft2(abs(fft2(M, nf(1), nf(2))).^2)), 1);
lag = -AHW:AHW;
R = acf(mod(lag, nf(1)) + 1, mod(lag, nf(2)) + 1); R = R / R(AHW+1, AHW+1);
[X, Y] = meshgrid(lag);
R = R .* cos(pi/2 * min(hypot(X, Y)/(AHW + 1), 1)).^2; %taper
Rp = zeros(nf); Rp(mod(lag, nf(1)) + 1, mod(lag, nf(2)) + 1) = R;
S = real(fft2(Rp)); S = max(S, SFLOOR*max(S(:))); %noise power spectrum
kl = -KHW:KHW;
Kp = zeros(nf); Kp(mod(kl, nf(1)) + 1, mod(kl, nf(2)) + 1) = K;
flt = struct('KS', conj(fft2(Kp)) ./ S, 'nf', nf, 'med', 0, 'sc', 1);
F = applyFilt(flt, Z);
x = F(~isnan(F)); flt.med = median(x); flt.sc = mad(x, 1)/0.6741891400433162;
F = (F - flt.med)/flt.sc;
end


function K = sourceKernel(Z, hw)
% average of patches around strong maxima of Z with no other maximum above 5
% within 6 px, each centred on its centroid; symmetric under a half turn (so
% rows and columns keep their own extents), cut to a disk of radius hw and
% normalized to unit sum. Falls back to weaker maxima, then to a Gaussian with
% the widths of the strong maxima along rows and columns, when there are too few
[H, W] = size(Z); [X, Y] = meshgrid(-hw:hw);
E = Z; E(isnan(E)) = -inf;
[r, c] = find(E == ordfilt2(E, 9, ones(3)) & E > 5);
v = E(sub2ind([H W], r, c));
D = hypot(r - r', c - c'); D(1:numel(r)+1:end) = inf;
inside = r > hw + 1 & r < H - hw - 1 & c > hw + 1 & c < W - hw - 1;
K = [];
for zMin = [15 10 8]
    k = find(v > zMin & min(D, [], 2) > 6 & inside);
    [~, o] = sort(v(k), 'descend'); k = k(o(1:min(end, 60)));
    acc = zeros(2*hw + 1); n = 0;
    for q = k'
        rr = r(q)-1:r(q)+1; cc = c(q)-1:c(q)+1;
        w = Z(rr, cc); w = max(w - min(w(:)), 0); w(isnan(w)) = 0;
        [CC, RR] = meshgrid(cc, rr);
        cen = [sum(RR(:).*w(:)), sum(CC(:).*w(:))] / sum(w(:));
        pt = interp2(Z, X + cen(2), Y + cen(1), 'linear');
        if any(isnan(pt(:))), continue, end
        acc = acc + pt/pt(hw+1, hw+1); n = n + 1;
    end
    if n >= 5
        K = acc/n;
        break
    end
end
if isempty(K)
    [sy, sx] = peakWidths(Z);
    K = exp(-(Y.^2/(2*sy^2) + X.^2/(2*sx^2)));
end
K = K - median([K(1, :), K(end, :), K(:, 1)', K(:, end)']);
K = (K + rot90(K, 2))/2;
K = K .* (hypot(X, Y) <= hw + 0.5);
K = K / sum(K(:));
end


function [sy, sx] = peakWidths(Z)
% Gaussian widths along rows and columns from the curvature of log Z at maxima
% above 10 (3 pixels per axis; median over at least 5 maxima, else 1 px)
[H, W] = size(Z); E = Z; E(isnan(E)) = -inf;
[r, c] = find(E == ordfilt2(E, 9, ones(3)) & E > 10);
in = r > 1 & r < H & c > 1 & c < W; r = r(in); c = c(in);
s = nan(numel(r), 2);
for k = 1:numel(r)
    p = [Z(r(k), c(k)), Z(r(k) - 1, c(k)), Z(r(k) + 1, c(k)), Z(r(k), c(k) - 1), Z(r(k), c(k) + 1)];
    if ~all(p > 0), continue, end
    ay = log(p(2)) + log(p(3)) - 2*log(p(1)); ax = log(p(4)) + log(p(5)) - 2*log(p(1));
    if ay < 0 && ax < 0, s(k, :) = [sqrt(-1/ay), sqrt(-1/ax)]; end
end
s = s(all(~isnan(s), 2), :);
sy = 1; sx = 1;
if size(s, 1) >= 5
    sy = min(max(median(s(:, 1)), 0.4), 3); sx = min(max(median(s(:, 2)), 0.4), 3);
end
end


function Z = combineTrials(zT, nanThresh)
% mean of the trials' z images, robust z; NaN where too few trials are imaged
Z = mean(zT, 3, 'omitnan');
Z(mean(isnan(zT), 3) > nanThresh) = nan;
Z = robustZ(Z);
end


function z = robustZ(a)
% z-score with median and MAD over the image, as getActImPeaks sets its threshold
v = a(~isnan(a));
z = (a - median(v)) / (mad(v, 1) / 0.6741891400433162);
end


function F = filtNan(A, K)
% normalized convolution ignoring NaNs; NaN stays NaN
m = ~isnan(A); A(~m) = 0;
F = imfilter(A, K, 0) ./ imfilter(double(m), K, 0);
F(~m) = nan;
end
