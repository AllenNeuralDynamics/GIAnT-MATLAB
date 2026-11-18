function [warp, info] = alignTrialsPatchwise(M, tforms, Mhalves, patchPx, nRounds)
%ALIGNTRIALSPATCHWISE Align trial mean images with a patchwise displacement field.
%
%   [WARP, INFO] = ALIGNTRIALSPATCHWISE(M, TFORMS0, MHALVES) aligns the trial
%   mean images M (H x W x T, NaN where unscanned) to a common frame, starting
%   from TFORMS0 (2 x 3 x T, e.g. rigid shifts, mapping the common frame
%   [row; col; 1] to each trial's image; trials whose transform is NaN are
%   skipped). A common-frame pixel x is read from trial t at x + u_t(x), with
%   u_t a smooth field set by its values d_kt at the centres c_k of patches of
%   imaged pixels:
%
%       u_t(x) = p_t(x) + sum_k phi_k(x) * (d_kt - p_t(c_k))
%
%   p_t is the plane (affine) least-squares fit through all the d_kt, and the
%   phi_k are Gaussians of SD WARP.sigma (half the median spacing of the
%   centres) around the centres, normalized to sum to 1. An affine field is
%   reproduced exactly; a centre's departure from the plane bends the field
%   around that centre. See patchFieldMatrix, trialWarpCoords.
%
%   Patches are k-means clusters of the coordinates of the pixels imaged in at
%   least half of the trials (under TFORMS0), about PATCHPX (default 1500)
%   imaged pixels each, so every centre is the centre of mass of a compact set
%   of imaged pixels, wherever those lie. Clusters with fewer than PATCHPX/2
%   pixels are dissolved into their neighbours. Patches much smaller than the
%   default can be dominated by a single spine, letting the field stretch or
%   shrink that spine.
%
%   The fit is groupwise, as in alignTrialsAffine: the template is the mean of
%   the warped trials, all the d_kt of a trial are refit to it jointly by
%   Gauss-Newton (a gain and offset per patch, sqrt images, NaNs ignored,
%   robust bisquare weights so that a spine changing brightness does not pull
%   the field), and the template is rebuilt, for NROUNDS rounds (default 3).
%   The mean field across trials is held at that of TFORMS0, so the common
%   frame keeps the trials' average geometry. A weak penalty puts centres with
%   little imaged structure in a trial on the plane.
%
%   MHALVES (H x W x T x 2) holds two independent halves of the data in M, of
%   equal size, so that their mean is M (e.g. interleaved blocks of frames,
%   balanced as loadAndProcessTrialAsync does). Each trial is refit on both halves;
%   their disagreement gives the noise of each centre's departure from the
%   plane, and the departures are scaled by max(0, 1 - noise variance /
%   variance across trials), so a patch whose bend is not above its noise
%   follows the plane. With MHALVES empty, the departures are kept as fit.
%
%   WARP has fields nodes (K x 2 centres, [row col]), sigma, disp (2 x K x T,
%   the d_kt as [row; col] in pixels), patch (H x W patch numbers, 0 outside),
%   and affine (2 x 3 x T, the plane p_t as a map like TFORMS0). If too few
%   pixels are imaged, WARP has no centres and affine is TFORMS0. INFO has nPx
%   (1 x K), noiseSD and bendSD (per axis, of each centre's departure from the
%   plane before scaling), weight (1 x K), corr (1 x T, each trial's
%   correlation with the final template over the patches), and rmsDisp (1 x T,
%   the RMS change of the field from TFORMS0 over the patches, in pixels).

if nargin < 4 || isempty(patchPx)
    patchPx = 1500;
end
if nargin < 5 || isempty(nRounds)
    nRounds = 3;
end
lambdaRel = 1e-2; %penalty on departures from the plane, relative to a typical centre's data term
[H, W, nT] = size(M, [1 2 3]);
trials = find(squeeze(all(~isnan(tforms), [1 2])))';
S = sqrtNan(M);
[cc, rr] = meshgrid(1:W, 1:H);
warp = struct('nodes', zeros(0,2), 'sigma', nan, 'disp', zeros(2,0,nT), 'patch', zeros(H,W), 'affine', tforms);
info = struct('nPx', [], 'noiseSD', [], 'bendSD', [], 'weight', [], 'corr', nan(1,nT), 'rmsDisp', nan(1,nT));
if numel(trials) < 2
    return
end

%patches, on the pixels imaged in at least half of the trials under TFORMS0
seen = zeros(H, W);
for t = trials
    seen = seen + ~isnan(warpImageAffine(S(:,:,t), tforms(:,:,t), 'linear'));
end
valid = seen >= numel(trials)/2;
if nnz(valid) < 500
    return
end
[patch, nodes] = makePatches(valid, patchPx);
K = size(nodes, 1);
if K > 1
    dd = sqrt((nodes(:,1) - nodes(:,1)').^2 + (nodes(:,2) - nodes(:,2)').^2);
    dd(1:K+1:end) = inf;
    sigma = median(min(dd, [], 2)) / 2;
else
    sigma = inf;
end
px = find(patch > 0);
r = rr(px); c = cc(px); pIx = patch(px);
[L, E] = patchFieldMatrix(nodes, sigma, r, c);

%initial centre displacements from TFORMS0
D = nan(2, K, nT);
for t = trials
    D(:,:,t) = tforms(:,1:2,t)*nodes' + tforms(:,3,t) - nodes';
end
D0 = D;
Dmean0 = mean(D(:,:,trials), 3);

for rnd = 1:nRounds
    template = buildTemplate(S, D, L, r, c, trials);
    for t = trials
        D(:,:,t) = refineField(S(:,:,t), template, D(:,:,t), L, r, c, pIx, E, lambdaRel);
    end
    D(:,:,trials) = D(:,:,trials) - mean(D(:,:,trials), 3) + Dmean0;
end

%reliability: refit each half against the final template
bend = nan(2, K, nT);
for t = trials
    bend(:,:,t) = (E*D(:,:,t)')';
end
bendVar = mean(reshape(permute(bend(:,:,trials), [2 1 3]), K, []).^2, 2)';
if isempty(Mhalves)
    noiseVar = nan(1, K);
    weight = ones(1, K);
else
    template = buildTemplate(S, D, L, r, c, trials);
    Dh = nan(2, K, nT, 2);
    for h = 1:2
        Sh = sqrtNan(Mhalves(:,:,:,h));
        for t = trials
            Dh(:,:,t,h) = refineField(Sh(:,:,t), template, D(:,:,t), L, r, c, pIx, E, lambdaRel);
        end
        Dh(:,:,trials,h) = Dh(:,:,trials,h) - mean(Dh(:,:,trials,h), 3) + Dmean0;
    end
    dBend = nan(2, K, nT);
    for t = trials
        dBend(:,:,t) = (E*(Dh(:,:,t,1) - Dh(:,:,t,2))')';
    end
    noiseVar = mean(reshape(permute(dBend(:,:,trials), [2 1 3]), K, []).^2, 2)' / 4; %halves of equal size (see MHALVES) each have twice the noise variance of the full data
    weight = max(0, 1 - noiseVar ./ bendVar);
    weight(~(bendVar > 0)) = 0;
end
for t = trials
    D(:,:,t) = D(:,:,t) - (1 - weight) .* bend(:,:,t);
end

warp.nodes = nodes;
warp.sigma = sigma;
warp.disp = D;
warp.patch = patch;
[~, ~, Lp] = patchFieldMatrix(nodes, sigma, [0; 1; 0], [0; 0; 1]); %plane at the origin and one pixel along each axis
for t = trials
    p = Lp*D(:,:,t)';
    warp.affine(:,:,t) = [eye(2) + [p(2,:) - p(1,:); p(3,:) - p(1,:)]', p(1,:)'];
end
info.nPx = accumarray(pIx, 1, [K 1])';
info.noiseSD = sqrt(noiseVar);
info.bendSD = sqrt(bendVar);
info.weight = weight;
template = buildTemplate(S, D, L, r, c, trials);
for t = trials
    u = L*D(:,:,t)';
    v = interp2(S(:,:,t), c + u(:,2), r + u(:,1), 'linear');
    ok = ~isnan(v) & ~isnan(template);
    info.corr(t) = corr(v(ok), template(ok));
    du = u - L*D0(:,:,t)';
    info.rmsDisp(t) = sqrt(mean(sum(du.^2, 2)));
end
end


function S = sqrtNan(M)
S = sqrt(max(M, 0)); %roughly variance-stabilized for shot noise
S(isnan(M)) = nan;
end


function [patch, nodes] = makePatches(valid, patchPx)
%k-means on the coordinates of the valid pixels, about patchPx pixels per
%cluster; while a cluster has fewer than patchPx/2 pixels, the smallest is
%dropped and k-means rerun from the remaining centres
[r, c] = find(valid);
X = [r c];
K = max(1, round(size(X,1) / patchPx));
opts = statset('MaxIter', 500, 'Streams', RandStream('mt19937ar', 'Seed', 0));
[lab, nodes] = kmeans(X, K, 'Replicates', 3, 'Options', opts);
n = accumarray(lab, 1, [K 1]);
while K > 1 && min(n) < patchPx/2
    [~, k] = min(n);
    nodes(k,:) = [];
    K = K - 1;
    [lab, nodes] = kmeans(X, K, 'Start', nodes, 'Options', statset('MaxIter', 500));
    n = accumarray(lab, 1, [K 1]);
end
patch = zeros(size(valid));
patch(valid) = lab;
end


function template = buildTemplate(S, D, L, r, c, trials)
%mean of the trials warped into the common frame at the fit pixels; pixels
%seen by fewer than half of the trials are NaN
V = nan(numel(r), numel(trials));
for k = 1:numel(trials)
    u = L*D(:,:,trials(k))';
    V(:,k) = interp2(S(:,:,trials(k)), c + u(:,2), r + u(:,1), 'linear');
end
template = mean(V, 2, 'omitnan');
template(sum(~isnan(V), 2) < numel(trials)/2) = nan;
end


function d = refineField(S, T, d, L, r, c, pIx, E, lambdaRel)
%robust Gauss-Newton for all centre displacements d (2 x K) of one trial:
%minimize the bisquare loss of S(x + u(x)) - a_k*T(x) - b_k over the fit
%pixels x (k their patch), with u = L*d', plus lambda*|E*d'|^2
maxIter = 30;
tol = 1e-3; %pixels, largest centre update
K = size(L, 2);
n = numel(r);
[gc, gr] = gradient(S);
Mp = sparse((1:n)', pIx, 1, n, K);
EE = E'*E;
w = ones(n, 1);
for it = 1:maxIter
    u = L*d';
    rq = r + u(:,1); cq = c + u(:,2);
    Wt = interp2(S, cq, rq, 'cubic');
    GR = interp2(gr, cq, rq, 'linear');
    GC = interp2(gc, cq, rq, 'linear');
    ok = ~isnan(T) & ~isnan(Wt) & ~isnan(GR) & ~isnan(GC);
    if nnz(ok) < 200
        return
    end
    Tz = T; Tz(~ok) = 0;
    y = Wt; y(~ok) = 0;
    GR(~ok) = 0; GC(~ok) = 0;
    wk = w .* ok;
    %Wt(d + dd) ~ Wt + GR.*(L*ddr) + GC.*(L*ddc) = a_k*T + b_k
    Xs = [spdiags(Tz, 0, n, n)*Mp, Mp];   %gain, offset per patch
    Xd = [-GR .* L, -GC .* L];           %centre updates, rows then cols
    WXd = wk .* Xd;
    A11 = full(Xs' * (wk .* Xs));
    A12 = full(Xs' * WXd);
    A22 = Xd' * WXd;
    b1 = full(Xs' * (wk .* y));
    b2 = WXd' * y;
    lambda = lambdaRel * median(diag(A22));
    A22 = A22 + lambda * blkdiag(EE, EE);
    b2 = b2 - lambda * [EE*d(1,:)'; EE*d(2,:)'];
    A = [A11, A12; A12', A22];
    A = A + (1e-10 * mean(diag(A)) + eps) * eye(4*K); %patches with no usable pixels in this trial
    beta = A \ [b1; b2];
    e = y - Xs*beta(1:2*K) - Xd*beta(2*K+1:end);
    sc = 1.4826 * accumarray(pIx(ok), abs(e(ok)), [K 1], @median) + eps;
    uu = e(ok) ./ (4.685 * sc(pIx(ok)));
    w(ok) = (1 - uu.^2).^2 .* (abs(uu) < 1);
    dd = [beta(2*K+1:3*K)'; beta(3*K+1:4*K)'];
    step = max(sqrt(sum(dd.^2, 1)));
    if step > 1
        dd = dd/step;
    end
    d = d + dd;
    if step < tol
        break
    end
end
end
