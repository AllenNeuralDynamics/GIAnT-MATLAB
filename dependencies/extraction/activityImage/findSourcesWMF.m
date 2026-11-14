function [seeds, info] = findSourcesWMF(actT, params)
%FINDSOURCESWMF Source seeds from the trials' aligned activity images.
%
%   [SEEDS, INFO] = FINDSOURCESWMF(ACTT, PARAMS) finds sources in ACTT
%   (rows x cols x trials), the trials' activity images in the common frame
%   (NaN where not imaged). SEEDS is N x 2 [row col], 1-indexed and subpixel,
%   strongest first. Uses params.wmfPeakth, params.wmfDip and params.nanThresh.
%
%   1. Each trial's image minus its local 5x5 median, as a robust z-score
%      (median and MAD over the image), as each trial's own peak finding does.
%   2. The combined image Z is the mean of the trials' z images, again as a
%      robust z-score; pixels imaged in fewer than 1 - nanThresh of the trials
%      are left out.
%   3. F is the whitened matched filter of Z: the detection statistic for a
%      source of the empirical shape K in noise of the empirical
%      autocorrelation R. K is the average of patches around strong, isolated
%      maxima of Z. R is the autocorrelation of Z away from activity, which
%      carries the correlation of the noise along scan lines. F is Z
%      correlated with R^-1 K, again as a robust z-score.
%   4. Local maxima of F above params.wmfPeakth are candidates. From the
%      strongest down, a candidate within MERGEPX of a stronger kept one is
%      merged into it unless F dips by at least params.wmfDip between them, and
%      by params.wmfDip/sqrt(2) in each half of the trials (alternate trials,
%      filtered alike). Two sources are kept apart only where they are resolved
%      reproducibly, so one source is not split into seeds at its edges.
%   5. A seed's position is the centroid of Z over the 3x3 pixels around its
%      maximum.
%   6. A seed must show activity (its trial's z image smoothed with a 1 px
%      Gaussian, above ZPRES within 1 px) in at least MINTRIALS trials.
%
%   With fewer than 2*MINHALF imaged trials (e.g. one long recording), the
%   halves are too noisy to test and a source may be active in only some of
%   so few trials: two maxima need only the dip in F to stay apart, and a seed
%   need only be active in one trial. On synthetic images this keeps more
%   sources at 2-3 trials than requiring activity in 2 of them, with no extra
%   seeds and no split sources.
%
%   INFO: Z, F, K, R, and per seed the value of F (value) and the number of
%   trials it is active in (nTrials).

MERGEPX = 6; ZPRES = 4; MINTRIALS = 2; MINHALF = 2;
[H, W, T] = size(actT);
seeds = zeros(0, 2); info = struct('Z', [], 'F', [], 'K', [], 'R', [], 'value', zeros(0, 1), 'nTrials', zeros(0, 1));
if T == 0 || all(isnan(actT(:)))
    return
end

zT = nan(H, W, T); zS = zT;
g = fspecial('gaussian', 7, 1);
for t = 1:T
    a = double(actT(:, :, t));
    if all(isnan(a(:))), continue, end
    zT(:, :, t) = robustZ(a - nanmedfilt2(a, [5 5]));
    zS(:, :, t) = robustZ(filtNan(zT(:, :, t), g));
end
imaged = find(squeeze(any(~isnan(zT), [1 2])))';
useHalves = numel(imaged) >= 2*MINHALF;
minTrials = MINTRIALS;
if ~useHalves
    minTrials = 1;
end
Z = combineTrials(zT, params.nanThresh);
[F, K, R] = wmfImage(Z);
if useHalves
    Fh = {wmfImage(combineTrials(zT(:, :, imaged(1:2:end)), params.nanThresh)), ...
          wmfImage(combineTrials(zT(:, :, imaged(2:2:end)), params.nanThresh))};
end

E = F; E(isnan(E)) = -inf;
[r, c] = find(E == ordfilt2(E, 9, ones(3)) & E > params.wmfPeakth);
v = E(sub2ind([H W], r, c));
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

ix = sub2ind([H W], round(seeds(:, 1)), round(seeds(:, 2)));
nTrials = zeros(numel(r), 1);
for t = 1:T
    z = zS(:, :, t); z(isnan(z)) = -inf;
    zm = ordfilt2(z, 9, ones(3));
    nTrials = nTrials + (zm(ix) > ZPRES);
end
ok = nTrials >= minTrials;
seeds = seeds(ok, :);
info = struct('Z', Z, 'F', F, 'K', K, 'R', R, 'value', v(ok), 'nTrials', nTrials(ok));
end


function [F, K, R] = wmfImage(Z)
% whitened matched filter of the combined z image Z (NaN where not used)
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
Z0 = Z; Z0(~m) = 0;
F = real(ifft2(fft2(Z0, nf(1), nf(2)) .* conj(fft2(Kp)) ./ S));
F = F(1:H, 1:W); F(~m) = nan;
F = robustZ(F);
end


function K = sourceKernel(Z, hw)
% average of patches around strong maxima of Z with no other maximum above 5
% within 6 px, each centred on its centroid; symmetrized, cut to a disk of
% radius hw and normalized to unit sum. Falls back to weaker maxima, then to a
% Gaussian (sigma 1 px), when there are too few
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
    K = exp(-(X.^2 + Y.^2)/2);
end
K = K - median([K(1, :), K(end, :), K(:, 1)', K(:, end)']);
K = (K + rot90(K, 2))/2;
K = K .* (hypot(X, Y) <= hw + 0.5);
K = K / sum(K(:));
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
