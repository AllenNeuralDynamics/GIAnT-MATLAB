function [tforms, info] = alignTrialsAffine(M, tforms, nRounds)
%ALIGNTRIALSAFFINE Refine per-trial alignment of trial mean images to an affine model.
%
%   [TFORMS, INFO] = ALIGNTRIALSAFFINE(M, TFORMS0) refines the transforms
%   TFORMS0 (2 x 3 x T, e.g. rigid shifts) that map the common frame to each
%   trial's image: [row_trial; col_trial] = TFORMS(:,:,t) * [row; col; 1].
%   M is H x W x T, one mean image per trial (NaN where unscanned). Trials
%   whose transform is NaN are skipped.
%
%   Within-trial registration is rigid, but trial mean images differ by a
%   smooth deformation that tracks motion (stretches of ~0.5 px per 100 px).
%   This fits it groupwise: the template is the mean of the trials warped
%   into the common frame, each trial's affine is refit to it by Gauss-Newton
%   (with a gain and offset, on sqrt-transformed images, ignoring NaNs), and
%   the template is rebuilt, for NROUNDS rounds (default 3). The mean linear
%   part across trials is removed each round, so the common frame keeps the
%   trials' average geometry, with the image centre where it was.
%
%   INFO.corr is each trial's correlation with the final template, and
%   INFO.rmsDisp the RMS displacement of the affine refinement relative to
%   TFORMS0 over the template's valid pixels, in pixels.

if nargin < 3 || isempty(nRounds)
    nRounds = 3;
end
[H, W, nT] = size(M, [1 2 3]);
trials = find(squeeze(all(~isnan(tforms), [1 2])))';
tforms0 = tforms;
S = sqrt(max(M, 0)); %roughly variance-stabilized for shot noise
S(isnan(M)) = nan; %max(NaN,0) is 0; unscanned pixels must stay NaN, not become sharp zero edges
[cc, rr] = meshgrid(1:W, 1:H);
x0 = [(H+1)/2; (W+1)/2]; %image centre, kept fixed when re-centring
info.corr = nan(1, nT);
info.rmsDisp = nan(1, nT);
if numel(trials) < 2
    return
end

for rnd = 1:nRounds
    template = buildTemplate(S, tforms, trials, rr, cc);
    for t = trials
        tforms(:,:,t) = refineAffine(S(:,:,t), template, tforms(:,:,t), rr, cc);
    end
    tforms = recentre(tforms, trials, x0);
end

template = buildTemplate(S, tforms, trials, rr, cc);
sel = ~isnan(template);
for t = trials
    w = warpImageAffine(S(:,:,t), tforms(:,:,t), 'linear');
    v = sel & ~isnan(w);
    info.corr(t) = corr(w(v), template(v));
    d = tforms(:,:,t) - tforms0(:,:,t);
    dr = d(1,1)*rr(sel) + d(1,2)*cc(sel) + d(1,3);
    dc = d(2,1)*rr(sel) + d(2,2)*cc(sel) + d(2,3);
    info.rmsDisp(t) = sqrt(mean(dr.^2 + dc.^2));
end
end


function template = buildTemplate(S, tforms, trials, rr, cc)
%mean of the trials warped into the common frame; pixels seen by fewer than
%half of the trials are NaN
stack = nan(size(S,1), size(S,2), numel(trials));
for k = 1:numel(trials)
    stack(:,:,k) = warpImageAffine(S(:,:,trials(k)), tforms(:,:,trials(k)), 'linear');
end
template = mean(stack, 3, 'omitnan');
template(sum(~isnan(stack), 3) < numel(trials)/2) = nan;
end


function A = refineAffine(S, template, A, rr, cc)
%forward-additive Gauss-Newton: minimize sum over valid pixels of
%(S(A*x) - a*template(x) - b)^2 over A, a, b
maxIter = 50;
tol = 1e-3;  %pixels, largest displacement change over the template
sel = find(~isnan(template));
if numel(sel) < 500
    return
end
r = rr(sel); c = cc(sel); T = template(sel);
r0 = mean(r); c0 = mean(c); s = 100; %centred, scaled coordinates for conditioning
u = (r - r0)/s; v = (c - c0)/s;
[gc, gr] = gradient(S);
for it = 1:maxIter
    rT = A(1,1)*r + A(1,2)*c + A(1,3);
    cT = A(2,1)*r + A(2,2)*c + A(2,3);
    Wt = interp2(S, cT, rT, 'cubic');
    GR = interp2(gr, cT, rT, 'linear');
    GC = interp2(gc, cT, rT, 'linear');
    ok = ~isnan(Wt) & ~isnan(GR) & ~isnan(GC);
    if nnz(ok) < 500
        return
    end
    %Wt(A + dA) ~ Wt + J*dp, with dp the update in centred, scaled coordinates
    J = [GR(ok).*u(ok), GR(ok).*v(ok), GR(ok), GC(ok).*u(ok), GC(ok).*v(ok), GC(ok)];
    beta = [T(ok), ones(nnz(ok),1), -J] \ Wt(ok);
    dp = beta(3:8);
    B = [dp(1) dp(2); dp(4) dp(5)];       %linear part, per scaled unit
    dL = B/s;
    dt = [dp(3); dp(6)] - dL*[r0; c0];
    %limit the step to 1 px anywhere on the template
    corners = [min(r) max(r) min(r) max(r); min(c) min(c) max(c) max(c)];
    stepMax = max(sqrt(sum((dL*corners + dt).^2, 1)));
    if stepMax > 1
        dL = dL/stepMax; dt = dt/stepMax;
    end
    A = A + [dL dt];
    if stepMax < tol
        break
    end
end
end


function tforms = recentre(tforms, trials, x0)
%right-multiply every transform by B so that the trials' mean linear part is
%the identity and the mean mapping of the image centre x0 is unchanged
Abar = mean(tforms(:,:,trials), 3);
Lbar = Abar(:,1:2);
Linv = inv(Lbar);
B = [Linv, x0 - Linv*x0; 0 0 1];
for t = trials
    Ah = [tforms(:,:,t); 0 0 1] * B;
    tforms(:,:,t) = Ah(1:2,:);
end
end
