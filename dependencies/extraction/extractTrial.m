function varargout = extractTrial(Y_obs,Finv, sources, selPix, params, alignData, GT)

Y_obs = double(Y_obs);
Finv = double(Finv);
sz = size(selPix);

if isempty(params.photonScale) %If not provided, estimate the standard deviation of a dim pixel
    pxSTD = nan(1,size(Y_obs,1));
    mY = nan(1,size(Y_obs,1));
    for pxIx = 1:size(Y_obs,1)
        [pxSTD(pxIx), mY(pxIx)] = std(Y_obs(pxIx,:,1),1./Finv(pxIx,:),2,'omitmissing');
    end
    Vb = prctile(pxSTD.^2, 5,'all');
    selBright = mY(:) > prctile(mY(:), 40) & mY(:) < prctile(mY(:), 90);
    params.photonScale = prctile((pxSTD(selBright).^2-Vb)./mY(selBright), 10);
end

%rescale data
Y_obs = Y_obs./params.photonScale;
if nargin>6
    GT.B = GT.B./params.photonScale;
    GT.S = GT.S./params.photonScale;
    GT.X = GT.X./params.photonScale;
    GT.Y = GT.Y./params.photonScale;
end
% params.lambda = 1; %after the above normalization, setting this higher than 1 (e.g. 2) encourages sparsity of activity

%break the problem into separable chunks, i.e. nonoverlap sets
zones = false(sz); zones(sub2ind(sz,sources.R,sources.C)) = true;
zones = imdilate(zones, strel('disk', ceil(1.5*params.sigma_px+(size(params.validKernel,1)-1)/2)));
CC = bwconncomp(zones,4);
nProblems = CC.NumObjects;

selIdxs = nan(sz);
selIdxs(selPix) = 1:sum(selPix(:));

if nargout==1
    doParallel = true;
    if isempty(gcp('nocreate'))
        parpool('Threads'); % start pool if none
    end
    analysisFutures = parallel.FevalFuture.empty(nProblems,0);
else
    doParallel = false;
end

sourceList = cell(1,nProblems); pxList_p = cell(1,nProblems);
for problemIx = 1:nProblems %9
    disp(['Solving subproblem ' int2str(problemIx) ' of ' int2str(CC.NumObjects)]);
    pxList = CC.PixelIdxList{problemIx};
    
    selSources = any(sub2ind(sz,sources.R,sources.C)==pxList',2); 
    sources_p.R = sources.R(selSources);
    sources_p.C = sources.C(selSources);
    sourceList{problemIx} = find(selSources);
    
    selPix_p = false(sz);
    selPix_p(sub2ind(sz,sources_p.R,sources_p.C)) = true;
    selPix_p = imdilate(selPix_p, strel('disk',params.selRadius)) & selPix;
    pxList_p{problemIx} = find(selPix_p(:));

    Y_p = Y_obs(selIdxs(pxList_p{problemIx}),:,:); %observations
    F_inv_p = Finv(selIdxs(pxList_p{problemIx}),:); % inverse freshness

    %crop, to reduce memory and improve visualization
    rowSupport = [find(any(selPix_p,2),1,'first') find(any(selPix_p,2),1,'last')];
    colSupport = [find(any(selPix_p,1),1,'first') find(any(selPix_p,1),1,'last')];
    sources_p.R = sources_p.R-rowSupport(1)+1;
    sources_p.C = sources_p.C-colSupport(1)+1;
    selPix_p = selPix_p(rowSupport(1):rowSupport(2), colSupport(1):colSupport(2));

    if nargin>6 %Ground truth supplied
        GTp.S = GT.S(selSources,:); % spikes,sources x time
        GTp.X = GT.X(selSources,:); % S convolved with kernel
        GTp.B = GT.B(selIdxs(pxList_p{problemIx}),:); % background
        GTp.H = GT.H(selIdxs(pxList_p{problemIx}), selSources); % source images; pixels x sources
        GTp.Hs = GT.Hs(selIdxs(pxList_p{problemIx}), selSources); % superres source images; pixels x sources

        if doParallel
            analysisFutures(problemIx) = parfeval(@extractSources, 6, Y_p, F_inv_p, sources_p, selPix_p, params, GTp);
        else
            [Hi,Si,Bi, LSi, SNRi, errFinal] = extractSources(Y_p, F_inv_p, sources_p, selPix_p, params, GTp); %#ok<ASGLU>
        end
    else
        if doParallel
            analysisFutures(problemIx) = parfeval(@extractSources, 5, Y_p, F_inv_p, sources_p, selPix_p, params);
        else
            [Hi,Si,Bi, LSi, SNRi] = extractSources(Y_p, F_inv_p, sources_p, selPix_p, params); %#ok<ASGLU>
        end
    end

    if ~doParallel %process this problem now
        %TODO: assemble Hi/Si/Bi/LSi/SNRi into H/B/S/LS/F0/SNR across problems
        error('extractTrial:NonParallelNotImplemented', ...
            ['Non-parallel extractTrial assembly is not implemented. ', ...
             'Call with a single output argument to use the parallel futures path.']);
    end
end

if doParallel
    % Capture scale used for this trial (may have been estimated above).
    photonScale = params.photonScale;
    varargout{1} = afterAll(analysisFutures,@(x)assembleResults(x, pxList_p, selIdxs, sourceList, numel(sources.R), size(Y_obs,2), photonScale),6, "PassFuture",true);
else
    varargout = {H,B,S,LS,F0,SNR}; %#ok<USENS>
end
end

function [H,B,S,LS,F0,SNR] = assembleResults(analysisFutures, pxList_p, selIdxs, sourceList, nSources, nTimepoints, photonScale)
%assemble results
%sz = CC.ImageSize;
nFutures = length(analysisFutures);
%H = nan([sz nSources], 'single');
for j = 1:nFutures
    [idx, Hi,Si,Bi, LSi, SNRi] = fetchNext(analysisFutures);
    if j==1
        nChan = size(Bi,3);
        H = nan(nnz(selIdxs>0),nSources, 'single');
        B = nan(nnz(selIdxs>0), nTimepoints, nChan);
        S = nan(nSources, nTimepoints, nChan);
        LS = nan(nSources, nTimepoints, nChan);
        F0 = nan(nSources, nTimepoints, nChan);
        SNR = nan(nSources,1);
    end

    %nj = size(Hi,2); %number of sources in this subproblem
    pxList = pxList_p{idx};
    
    %tmp = nan([prod(sz), nj]);
    %tmp(pxList,:) = Hi;
    %H(:,:,sourceList{idx}) = reshape(tmp, [sz nj]); %footprints
    H(selIdxs(pxList),sourceList{idx}) = Hi;
    B(selIdxs(pxList),:,:) = Bi; %background
    
    S(sourceList{idx},:,:) = Si; %spikes
    LS(sourceList{idx},:,:) = LSi; %least squares estimate
    for chIx = 1:size(Bi,3)
        F0(sourceList{idx},:,chIx) = (Hi'*Bi(:,:,chIx))./sum(Hi.^2,1)'; %F0 is weighted background
    end
    SNR(sourceList{idx}) = SNRi;
end

% Restore movie fluorescence units (optimization ran in photon-normalized space).
% H and SNR are scale-invariant / sum-normalized and are left unchanged.
B = B .* photonScale;
S = S .* photonScale;
LS = LS .* photonScale;
F0 = F0 .* photonScale;
end




function [H_est,S_est_new, B_est, dFls, Xsnr, errFinal] = extractSources(Y_obs, Finv, sources, selPix, params, GT)
%performs source extraction by alternating coordinate optimization on the
%footprints (H,Hs), baseline (B), and source activities (S,X), in a constrained NMF framework

%Yobs is [#selected pixels]x[#timepoints]
%Finv is [#selected pixels]x[#timepoints], 1./freshness, "freshness" is proportional to #measurements averaged into each sample
%selpix is a 2D logical map, containing [#selected pixels] true values, that maps Yobs to 2D space
%params

%outputs the superresolution source estimate (Hs), the superresolution
%spikes (S), and the baseline (B)

num_channels = size(Y_obs,3);
if num_channels>1
    Y2 = Y_obs(:,:,2);
    Y_obs = Y_obs(:,:,1);
else
    Y2 = [];
end

sz = size(selPix);
num_sources = numel(sources.R);
num_time_points = size(Y_obs,2);

%Initial estimates for B, S, H
for six = num_sources:-1:1
    tmp = zeros(sz);
    tmp(sources.R(six), sources.C(six)) = 1;

    tmpValid = imdilate(logical(tmp), params.validKernel);

    tmpHs = imdilate(tmp, strel('disk', 1)); tmpHs = tmpHs./sum(tmpHs,'all');
    tmpH = imgaussfilt(tmpHs, params.sigma_px);
    
    Hs_est(:,six) = tmpHs(selPix); %initial estimate
    H_est(:,six) = tmpH(selPix); %initial estimate
    Hvalid(:,six) = tmpValid(selPix); %spatial support for each source in solver
end

%initialize B
if ~isfield(params, 'denoiseWindow_samps')
    params.denoiseWindow_samps = params.denoiseWindow_s .* params.analyzeHz;
end
params.denoiseWindow_samps = ceil(params.denoiseWindow_samps);
%denoised = smoothdata(Y_obs,2,"movmean",,'omitmissing');
B_est = max(params.minBaseline, splitFreq(Y_obs, params.denoiseWindow_samps, ceil(params.baselineWindow_samps/params.denoiseWindow_samps)));

%medRes = median(denoised-LP,2);
typicalX = sqrt(mean((Y_obs(:,1:100)-B_est(:,1:100)).^2,'all'))*ones(num_sources,num_time_points);

dFls = H_est\(Y_obs-B_est);

%initialize S
S_est = max(0, dFls - [zeros(size(dFls,1),1) dFls(:,1:end-1)])./params.k(floor(end/2)+1);

%overwrite with GT for testing?
% B_est = GT.B;
% S_est = GT.S;

problemS.lb = -eps*ones(size(S_est));
problemS.ub = inf(size(S_est));

problemH.lb = -eps*ones(size(Hs_est));
problemH.ub = eps*ones(size(Hs_est));
problemH.ub(Hvalid) = inf;

if nargin>6 %ground truth supplied
    B_err = nan(1, params.nmfIter);
    H_err = nan(1, params.nmfIter);
    S_err = nan(1, params.nmfIter);
    B_err(1) = mean(B_est(:)-GT.B(:)).^2;
    H_err(1) = mean(H_est(:)-GT.H(:)).^2;
    S_err(1) = 1-corr(S_est(:),GT.S(:));
    hAx = plotGT([], S_est, Hs_est, B_est, S_err, H_err, B_err, GT, selPix, params);
end

%optimization options
opts = optimoptions('fmincon', ...
    'Algorithm','trust-region-reflective', ...   %'interior-point'
    'SpecifyObjectiveGradient',true,'SubproblemAlgorithm', 'cg', ...
    'Display','none','MaxPCGIter', 10);

doFitS = true;
doFitH = true;
doFitB = true;

%OPTIMIZE
for outerLoop = 1:params.nmfIter
    opts.MaxIterations = 10*outerLoop;

    %SOLVE FOR S
    if doFitS
        opts.TypicalX = typicalX;
        % Objective function handle (returns [f,g, Hinfo])
        objS = @(x) objfun_S_wrapper(x, Y_obs, H_est, B_est, params.k, Finv, params.lambda);

        % Hessian multiply for fmincon signature (x,y,flag)
        opts.HessianMultiplyFcn = @(hinfo, v, flag) hessmult_S_wrapper(hinfo, Y_obs, H_est, B_est, params.k, Finv, params.lambda, v); %hessmult_S_wrapper takes (Svec, Z, H, B, k, F, lambda, v)

        % Call fmincon
        [S_est_new, lossS] = fmincon(objS, S_est, [], [], [], [], problemS.lb, problemS.ub, [], opts);
    else
        S_est_new = S_est;
    end

    if nargin>6 %if ground truth was supplied, make error plots
        B_err(outerLoop+1) = mean(B_est(:)-GT.B(:)).^2;
        H_err(outerLoop+1) = mean(H_est(:)-GT.H(:)).^2;
        S_err(outerLoop+1) = 1-corr(S_est_new(:),GT.S(:));
        plotGT(hAx, S_est_new, H_est, B_est, S_err, H_err, B_err, GT, selPix, params);
        errFinal = S_err(outerLoop+1);
    else
        errFinal = nan;
    end

    %update X
    X_est_new = convn(S_est_new, params.k, 'same');

    %estimate noise,SNR
    resid = Y_obs - (B_est + H_est*X_est_new);
    
    %compute weighted rms error estimate at each pixel
    residWeights = 1./Finv;
    residWeights(resid>=0) = 0;
    residVar = sum(resid.^2.*residWeights,2)./sum(residWeights,2);
    W = diag(1./residVar);
    covX = inv(H_est' * W * H_est); %uncertainty estimate for X
    Xnoise = sqrt(diag(covX)./params.tau_full);
    Xsnr = std(X_est_new, 0,2)./Xnoise;

    if outerLoop==params.nmfIter
        break %we stop optimizing after fitting S on last loop
    end

    Xfloor = computeFloor(X_est_new, params.denoiseWindow_samps, params.baselineWindow_samps);
    X_est_new = max(0, X_est_new - Xfloor - Xnoise);

    %SOLVE FOR B
    if doFitB
        B_est_new = fitB(B_est, Y_obs,X_est_new,H_est, Finv, selPix, params); %(Y,X, H, Finv, selPix, params)
        B_est = B_est_new;
    end

    %SOLVE FOR Hs
    if doFitH
        opts.TypicalX = max(H_est, 0.01);
        % Objective function handle (returns [f,g, Hinfo])
        objHs = @(Hs) objfun_Hs_wrapper(Hs, Y_obs, X_est_new, B_est, params.Hfilter, selPix, Finv, params.lambda); %(Hs_vec, Z, X, B, kk, selPix, F, lambda, v)
        opts.HessianMultiplyFcn = @(Hinfo, v, flag) hessmult_Hs_wrapper(Hinfo, Y_obs, X_est_new, B_est, params.Hfilter,selPix, Finv, params.lambda, v);
        % Call fmincon
        [Hs_est_new,lossH] = fmincon(objHs, Hs_est, [], [], [], [], problemH.lb, problemH.ub, [], opts);
    else
        Hs_est_new = Hs_est;
    end

    %update H
    tmp = zeros([sz num_sources]);
    tmp(repmat(selPix,1,1,num_sources)) = Hs_est_new;
    tmp = reshape(convn(tmp,params.Hfilter,'same'), numel(selPix), num_sources);
    H_est_new = tmp(selPix,:);

    %normalize new H and S
    normFac = sum(H_est_new, 1);

    H_est_new = H_est_new./normFac;
    Hs_est_new = Hs_est_new./normFac;
    S_est_new = S_est_new.*normFac';

    %update with new values
    H_est = H_est_new;
    Hs_est = Hs_est_new;
    S_est = S_est_new;
end

% Debias: freeze small events, refit free set (phi=0 => no L1). Avoids
% lb≈ub pinching, which breaks trust-region-reflective (trdog/quad1d).
free = S_est_new >= 1e-1;
S_est_new(~free) = 0;
if any(free(:))
    free_idx = find(free);
    szS = size(S_est_new);
    x0 = S_est_new(free);
    lambdaDebias = params.lambda * params.phi;
    optsDebias = opts;
    optsDebias.MaxIterations = 10*params.nmfIter;
    optsDebias.TypicalX = typicalX(free);
    optsDebias.HessianMultiplyFcn = @(hinfo, v, varargin) ...
        hessmult_S_free(hinfo, v, free_idx, szS, Y_obs, H_est, B_est, params.k, Finv, lambdaDebias);
    objS = @(x) objfun_S_free(x, free_idx, szS, Y_obs, H_est, B_est, params.k, Finv, lambdaDebias);
    x = fmincon(objS, x0, [], [], [], [], -eps*ones(size(x0)), inf(size(x0)), [], optsDebias);
    S_est_new(free) = x;
end

%update X
X_est_new = convn(S_est_new, params.k, 'same');

%estimate noise,SNR
resid = Y_obs - (B_est + H_est*X_est_new);

%compute weighted rms error estimate at each pixel
residWeights = 1./Finv;
residWeights(resid>=0) = 0;
residVar = sum(resid.^2.*residWeights,2)./sum(residWeights,2);
W = diag(1./residVar);
covX = inv(H_est' * W * H_est); %uncertainty estimate for X
Xnoise = sqrt(diag(covX)./params.tau_full);
Xsnr = std(X_est_new, 0,2)./Xnoise;

if any(imag(Xsnr(:)) ~= 0)
    warning('Xsnr is complex, taking real component');
    Xsnr = real(Xsnr);
end

%fit least squares
dFls = H_est\(Y_obs-B_est);

if ~isempty(Y2) % two-channel recording, process calcium data with same source footprints
    %fit initial baseline
    B2 = max(params.minBaseline, splitFreq(Y2, params.denoiseWindow_samps, ceil(params.baselineWindow_samps/params.denoiseWindow_samps)));
    opts.HessianMultiplyFcn = @(hinfo, v, flag) hessmult_S_wrapper(hinfo, Y2, H_est, B2, params.k2, Finv, params.lambda, v); %hessmult_S_wrapper takes (Svec, Z, H, B, k, F, lambda, v)
    opts.MaxIterations = 15;
    LS2 = H_est\(Y2-B2);
    objS = @(x) objfun_S_wrapper(x, Y2, H_est, B2, params.k2, Finv, params.lambda);
    [S2, ~] = fmincon(objS, LS2, [], [], [], [], problemS.lb, problemS.ub, [], opts);
    
    %X2 = convn(S2, params.k2, 'same'); %update X2
    
    % %estimate noise,SNR
    % resid = Y2- (B2 + H_est*X2);
    % residWeights = 1./Finv;
    % residWeights(resid>=0) = 0;
    % residVar = sum(resid.^2.*residWeights,2)./sum(residWeights,2);
    % W = diag(1./residVar);
    % covX = inv(H_est' * W * H_est); %uncertainty estimate for X
    % X2noise = sqrt(diag(covX)./params.tau2_samps);
    % X2snr = std(X2, 0,2)./X2noise;
    % X2floor = computeFloor(X2, params.denoiseWindow_samps, params.baselineWindow_samps);
    % X2 = max(0, X2 - X2floor - X2noise);
    % 
    % %Refit
    % B2 = fitB(B2, Y2,X2,H_est, Finv, selPix, params); %(Y,X, H, Finv, selPix, params)
    % LS2 = H_est\(Y2-B2);
    % objS = @(x) objfun_S_wrapper(x, Y2, H_est, B2, params.k2, Finv, params.lambda);
    % [S2, ~] = fmincon(objS, LS2, [], [], [], [], problemS.lb, problemS.ub, [], opts);
    % X2 = convn(S2, params.k2, 'same'); %update X2

    %append channel 2 S,B,LS, and SNR as additional dimension
    S_est_new = cat(3, S_est_new, S2);
    B_est = cat(3, B_est, B2); 
    dFls = cat(3, dFls, LS2); 
    %Xsnr = cat(2, Xsnr, X2snr);
end

end


function preds = sinePredictors(t,basePeriod)
T = 0:t-1;
maxN = ceil(t/basePeriod);
periods = (1:maxN)' .* basePeriod;
phase = 2*pi*T./periods;
preds = [sin(phase); cos(phase); T./max(T)-0.5];
end

function B = fitB(B0,Y,X, H, Finv, selPix, params)
%fits a baseline (i.e. F0) to the measurements in Y, given the current
%estimate of the activity, HX
num_time_points = size(Y,2);
[~, HP] = splitFreq(Y-H*X, 2*params.denoiseWindow_samps, ceil(params.baselineWindow_samps/params.denoiseWindow_samps));
HP(isnan(HP)) = 0;
HPsurround = getSurround(HP,selPix,params);
sinePreds = sinePredictors(num_time_points,params.baselineWindow_samps);
%sqrtF = sqrt(1./Finv);
%Yweighted = Y.*sqrtF; %multiply in the freshness effect
opts = optimoptions('lsqlin','Display','none');
for pxIx = size(Y,1):-1:1
    scale = max(1e-3, mean(Y(pxIx,:)));
    nValid = size(H,2);
    preds = [ H(pxIx,:)'.*X; sinePreds; squeeze(HPsurround(pxIx,:,:))'; ones(1,num_time_points).*scale;]';
    lb = -10*scale.*ones(1, size(preds,2)); lb(1:nValid) = 0;
    ub = 10*scale.*ones(1, size(preds,2));

    %calculate weighted regression
    W = sqrt(1./Finv(pxIx,:)) .*  B0(pxIx,:)./(B0(pxIx,:) + H(pxIx,:)*X);
    preds_scaled = W'.*preds;
    resp_scaled = W'.*Y(pxIx,:)';
    
    b = lsqlin(preds_scaled,resp_scaled,[],[],[],[],lb,ub,[],opts);
    B(pxIx,:) = max(params.minBaseline, (preds(:,(nValid+1):end)*b((nValid+1):end))');
end
end

function [LP,HP] = splitFreq(A, denoiseWindow, LPfactor)
nPages= floor(size(A,2)./denoiseWindow);
totSamps = nPages*denoiseWindow;
t = 1:size(A,2);

a = reshape(A(:,1:totSamps), size(A,1), denoiseWindow, nPages); % #pixels x #samps/page x #pages
Ma = squeeze(mean(a,2, 'omitmissing'));
SMa = smoothdata(Ma,2, 'lowess',LPfactor, 'omitmissing');
resid = Ma-SMa;
lowVals = resid<=ordfilt2(resid, max(2,ceil(0.15*LPfactor)), ones([1 LPfactor]));
Ma(~lowVals) = nan;
SMa = smoothdata(Ma,2, 'movmean',LPfactor, 'omitmissing');
for iter = 1:3
    selNans = isnan(SMa);
    if any(selNans(:))
        tmp = smoothdata(SMa,2, 'movmean',LPfactor, 'omitmissing');
        SMa(selNans) = tmp(selNans); 
    else
        break
    end
end
tDS = (denoiseWindow+1)/2 + (denoiseWindow).*(0:size(SMa,2)-1);
LP = nan(size(A));
for pxIx = 1:size(LP,1)
    LP(pxIx,:) = interp1(tDS,SMa(pxIx,:)',t,'linear','extrap');
end
HP = A-LP;
end

function HPsurround = getSurround(HP, selPix, params)
%convert to pixel space
HPfull = zeros(size(selPix,1), size(selPix,2), size(HP,2));
HPfull(repmat(selPix, 1, 1, size(HP,2))) = HP;
for filtIx = size(params.Bfilter,3):-1:1
    tmp = convn(HPfull, params.Bfilter(:,:,filtIx), 'same');
    tmp = reshape(tmp, [], size(tmp,3));
    HPsurround(:,:,filtIx) = tmp(selPix(:),:);
end
end

function hAx = plotGT(hAx, S_est, H_est, B_est, S_err, H_err, B_err, GT, selPix, params)
if isempty(hAx) % SET UP PLOTS
    figure;
    hAx(1) = subplot(3,2,1);
    hAx(2) = subplot(3,2,2); plot(B_err); title(hAx(2), 'B_err History');
    hAx(3) = subplot(3,2,3); title('S')
    hAx(4) = subplot(3,2,4); plot(S_err); title('S_err History')
    hAx(5) = subplot(3,2,5); title('H')
    hAx(6) = subplot(3,2,6); plot(H_err); title('H_err History')
end
sz = size(selPix);
cla(hAx(1))
imagesc(B_est - GT.B, 'parent', hAx(1), [-0.1 0.1]); title(hAx(1), '{B_{est} - B_{GT}}'); colorbar(hAx(1));
plot(hAx(2), B_err); title(hAx(2), 'B error history');
cla(hAx(3));
plot(conv2(GT.S, params.k, 'same')' + (1:size(GT.S,1)), 'r','parent', hAx(3)),
hold(hAx(3), 'on'),
plot(conv2(S_est, params.k, 'same')' + (1:size(GT.S,1)), 'b', 'parent', hAx(3));
set(hAx(3), 'xlim', [0 1000]);
title(hAx(3), '{X_{est} (blue) , X_{GT} (red)} (first 1000 samples)');
plot(hAx(4), S_err); title(hAx(4), 'S error history');
cla(hAx(5))
Him = nan(sz); Him(selPix) = sum(H_est,2);
imagesc(Him, 'parent', hAx(5)); title(hAx(5), 'H');
plot(hAx(6), H_err); title(hAx(6), 'H error history');
drawnow
end



function [f, gHs, HvHs] = objfun_Hs(Hs, Z, X, B, kk, selPix, F, lambda, v)
% OBJFUN_HS  loss, gradient, and Hessian-vector product w.r.t Hs
%
% Usage:
%   [f, gHs] = objfun_Hs_new(Hs, Z, X, B, kk, F, lambda);
%   [~, ~, HvHs] = objfun_Hs_new(Hs, Z, X, B, kk, F, lambda, VHs);
%
% Inputs:
%   Hs     : m-by-p matrix (optimization variable)
%   Z, B, F: m-by-n matrices
%   X      : p-by-n matrix
%   kk     : 2D convolution kernel
%   lambda : scalar (regularizer)
%   VHs    : optional m-by-p perturbation for Hessian-vector product
%
% Outputs:
%   f      : scalar loss
%   gHs    : m-by-p gradient (dL/dHs)
%   HvHs   : m-by-p Hessian-vector product (if VHs provided)
ns = size(Hs,2);
selPix3D = repmat(selPix, 1,1,ns);

% Forward: H from Hs
tmp = zeros([size(selPix) ns]);
tmp(repmat(selPix, 1, 1, ns)) = Hs;
tmp = convn(tmp, kk, 'same');
tmp = reshape(tmp,[numel(selPix), ns]);
H = tmp(selPix(:),:);

T = H * X + B;               % m-by-n
E = Z - T;                   % m-by-n
s = T + 1;% lambda;              % m-by-n

% Hessian-vector product (BUT NOT F and G)
if nargin >= 9 && ~isempty(v)
    for rix = size(v,2):-1:1
        % Forward perturbation
        tmp = zeros(size(selPix3D));
        tmp(selPix3D) = v(:,rix);
        tmp = convn(tmp, kk, 'same');
        V = reshape(tmp(selPix3D), size(Hs));

        % Perturbation in T
        dT = V * X;                     % m-by-n
        % coef = 2 .* (Z + lambda).^2 ./ ( F .* (s.^3) );   % m-by-n
        coef = 2 .* (Z + 1).^2 ./ ( F .* (s.^3) );   % m-by-n
        dp = coef .* dT;                % m-by-n

        % Map back to H-space
        HvH = dp * X.';                % m-by-p

        %Backprop to Hs in real space
        tmp = zeros(size(selPix3D));
        tmp(selPix3D) = HvH;
        tmp = convn(tmp, rot90(kk,2), 'same');
        HvHs(:,rix) = tmp(selPix3D);
    end
    f = [];
    gHs = [];
else %loss and gradient
    % Compute loss
    denom = F .* s;              % m-by-n
    f = sum( (E.^2) ./ denom , 'all' );

    % Gradient wrt T
    p = - (2 .* E .* s + E.^2) ./ ( F .* (s.^2) );  % m-by-n

    % Gradient wrt H
    gH = p * X.';                % m-by-p

    % % Chain rule: gradient wrt Hs
    tmp = zeros([size(selPix) ns]);
    tmp(selPix3D) = gH;
    tmp = convn(tmp, rot90(kk,2), 'same');
    gHs = reshape(tmp(selPix3D), size(Hs));

    HvHs = [];
end
end

function [f, gHs, Hinfo] = objfun_Hs_wrapper(Hs, Z, X, B, kk, selPix, F, lambda)
[f, gHs] = objfun_Hs(Hs, Z, X, B, kk, selPix, F, lambda);
Hinfo = Hs;
% if f<0 || ~isreal(gHs) || ~isreal(f)
%     keyboard
% end
end

function HvHs = hessmult_Hs_wrapper(Hs, Z, X, B, kk, selPix, F, lambda, v)
[~,~,HvHs] = objfun_Hs(Hs, Z, X, B, kk, selPix, F, lambda, v);
% if ~isreal(HvHs)
%     keyboard
% end
end

%%%%%
%%%%%

function [f, gS, HvS] = objfun_S(S, Z, H, B, k, F, lambda, v)
% OBJFUN_S  loss, gradient, Hessian-vector product w.r.t. S
%
% S : p-by-n
% Z,B,F : m-by-n
% H : m-by-p
% k : 1-D kernel (row or column) used in convn(S,k,'same')
% lambda : scalar regularizer, roughly the scale of 1 photon
% VS : optional p-by-n perturbation for Hv (same size as S)
%
% Outputs:
%  f   : scalar loss
%  gS  : p-by-n gradient dL/dS
%  HvS : p-by-n Hessian-vector product (if VS provided), else []
%
% The loss:
%   T = H * X + B,  X = convn(S, k, 'same')
%   f = sum( (Z - T).^2 ./ (F .* (T + lambda)), 'all' )

% Forward
X = convn(S, k, 'same');   % p-by-n
T = H * X + B;             % m-by-n
E = Z - T;                 % m-by-n
s = T + 1; %lambda;            % m-by-n

% Hessian-vector product branch
if nargin >= 8 && ~isempty(v)
    for rix = size(v,2):-1:1
        V = reshape(v(:,rix), size(S,2), size(S,1))';

        % Forward map of perturbation through conv: dX = convn(VS, k, 'same')
        dX = convn(V, k, 'same');    % p-by-n

        % Perturbation in T: dT = H * dX
        dT = H * dX;                  % m-by-n

        % Coefficient: coef = 2*(Z + lambda).^2 ./ (F .* s.^3)
        % coef = 2 .* (Z + lambda).^2 ./ ( F .* (s.^3) );  % m-by-n
        coef = 2 .* (Z + 1).^2 ./ ( F .* (s.^3) );  % m-by-n

        % dp = coef .* dT
        dp = coef .* dT;              % m-by-n

        % Map back: tmp = H.' * dp  (p-by-n)
        tmp = H.' * dp;              % p-by-n

        % Back through conv adjoint: HvS = convn(tmp, kflip, 'same')
        conv_adj = convn(tmp, flip(k), 'same');  % p-by-n
        HvS(:,rix) = conv_adj(:);
    end
    f = [];
    gS = [];
else
    % Objective value
    denom = F .* s;            % m-by-n
    f = sum( (E.^2) ./ denom, 'all' ) + lambda .* sum(S(:));

    % Gradient wrt T (elementwise)
    % p_elem = dL/dT = - (2 E s + E.^2) ./ (F .* s.^2)
    p_elem = - (2 .* E .* s + E.^2) ./ ( F .* (s.^2) );  % m-by-n

    % Gradient wrt X: gX = H.' * p_elem   (p-by-n)
    gX = H.' * p_elem;        % p-by-n

    % Gradient wrt S: adjoint of convn => convn(gX, flip(k), 'same')
    gS = convn(gX, flip(k), 'same') + lambda;  % p-by-n

    HvS = [];
end
end

function [f, gS, Hinfo] = objfun_S_wrapper(S, Z, H, B, k, F, lambda)
[f, gS] = objfun_S(S, Z, H, B, k, F, lambda); %objfun_S takes (S, Z, H, B, k, F, lambda, v)
Hinfo = S; %Hessian depends on S
% if f<0 || ~isreal(gS) || ~isreal(f)
%     keyboard
% end
end

function HvS = hessmult_S_wrapper(S, Z, H, B, k, F, lambda, v)
% Svec : current S flattened
% p = size(H,2);
% n = size(Z,2);
% S = reshape(S, p, n); %perhaps unnecessary?

[~,~,HvS] = objfun_S(S, Z, H, B, k, F, lambda, v);
% if ~isreal(HvS)
%     keyboard
% end
end

function [f, g, Hinfo] = objfun_S_free(x, free_idx, szS, Z, H, B, k, F, lambda)
%OBJFUN_S_FREE Objective on free S entries only (scatter/gather into full S).
S = zeros(szS);
S(free_idx) = x;
[f, gS] = objfun_S(S, Z, H, B, k, F, lambda);
g = gS(free_idx);
Hinfo = S;
end

function Hv = hessmult_S_free(Hinfo, v, free_idx, szS, Z, H, B, k, F, lambda)
%HESSMULT_S_FREE Hessian-vector product on free S entries only.
nv = size(v, 2);
Vfull = zeros(prod(szS), nv);
for j = 1:nv
    tmp = zeros(szS);
    tmp(free_idx) = v(:, j);
    Vfull(:, j) = tmp(:);
end
HvFull = hessmult_S_wrapper(Hinfo, Z, H, B, k, F, lambda, Vfull);
Hv = zeros(numel(free_idx), nv);
for j = 1:nv
    tmp = reshape(HvFull(:, j), szS);
    Hv(:, j) = tmp(free_idx);
end
end

function Xfloor = computeFloor(X, denoiseWindow, baseline)
ord = ceil(0.1*baseline); % a percentile filter to remove overfit small spikes during iterations
Xmed = medfilt2(X, [1 2*ceil(denoiseWindow)+1],"symmetric");
Xmed_min = ordfilt2(Xmed, ord, ones(1,ceil(baseline)), 'symmetric');
Xfloor = smoothdata(Xmed_min, 2,"movmean",ceil(baseline),"omitmissing");
end