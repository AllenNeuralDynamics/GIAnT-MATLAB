function [skIm, P, acc] = localizeSources_vIM(IM, vIM, params, doPlot)
%inputs:
%IM:        3D recording, X x Y x Time
%aData:     alignment metadata
%outputs:
%skIm:      activity image (sum of squared local-maximum scores / (reg + observations))
%P:         peaks of the activity image
%acc:       the parts of skIm (num, count, reg) and the number of frames
%           scored (nFrames), for accumulating across trials
nTimePoints = size(IM,3);
tau = params.tau_s.*params.alignHz; %time constant in frames
params.tau_frames = tau;
sigma = params.sigma_px; %space constant in pixels
baselineWindow = ceil(params.baselineWindow_Glu_s.*params.alignHz);
denoiseWindow = ceil(params.denoiseWindow_s.*params.alignHz);
nans = isnan(IM);

sz = size(IM);
valid = mean(nans,3)<params.nanThresh; %a pixel must be imaged at least (1-nanThresh) of the time to be included
if ~any(valid)
    warning('Recording had no valid pixels; likely too much motion')
    P = [];
    skIm = nan(sz(1:2));
    acc = struct('num', zeros(sz(1:2)), 'count', zeros(sz(1:2)), 'reg', 300, 'nFrames', 0);
    return
end
if nargin<4
    doPlot = false;
end

if isempty(vIM) %for raster imaging
    vIM = ones(sz);
end

%initialize filtered image
IMf = IM; clear IM;
IMf(repmat(~valid, 1, 1, nTimePoints)) = nan;
nans = isnan(IMf);
nanFrac = mean(nans,3);
vIM(nans) = nan; %1000*mean(vIM(:,:, 1:min(end,400)), 'all', 'omitnan');

% %fill in missing values
% IMf= reshape(IMf, sz(1)*sz(2), []);
% incomplete = nanFrac>0 & nanFrac<1;
% IMs = nan(size(IMf));
% IMs(incomplete(:),:) = smoothdata(IMf(incomplete(:),:), 2, 'movmean', baselineWindow, 'omitnan');
% IMf(reshape(nans, size(IMf))) = IMs(reshape(nans, size(IMf))); clear IMs
% IMf = reshape(IMf, sz(1),sz(2), []);

%smooth the data at a timescale on which fluctuations look more
%gaussian, for computing variances
IMs = smoothdata(IMf./vIM, 3, 'movmean', ceil(denoiseWindow/2), 'omitnan');
vIM = 1./smoothdata(1./vIM, 3, 'movmean', ceil(denoiseWindow/2), 'omitnan');
IMs = IMs.*vIM;

%baseline estimate
IMb = smoothdata(IMs, 3, 'movmedian', baselineWindow, 'omitnan');

%estimate Vb and Vk, parameters for estimating variance from baseline brightness
% Vb: the variance of a 'dim' pixel due to electronic and dark noise
% Vk: the slope of the variance-brightness relationship
firstValidFrames = find(any(~nans, [1 2]),500, 'first');
varIM = var(IMs(:,:,firstValidFrames),0,3,"omitmissing");
varIM(nanFrac>0.4) = nan;
Vb = params.VIF*prctile(varIM, 10, 'all');
varPred = mean(IMb(:,:,firstValidFrames),3,'omitmissing').* mean(vIM(:,:,firstValidFrames),3,'omitmissing');
selBright = varPred>prctile(varPred(:), 90);
Vk = prctile((varIM(selBright)-(Vb/params.VIF))./varPred(selBright), 10);

%Highpass filter in time; This must occur before DoG to avoid edge artifacts
IMf = IMf - IMb; 

stdIM = sqrt(max(0,Vk.*IMb.*vIM)+Vb); %compute standard deviation
%divide by uncertainty to get a Z-score
IMf = IMf./stdIM;
clear IMb vIM

%time matched filter
gamma = exp(-1/tau);
mem = max(0,gamma*IMf(:,:,end));
for t = size(IMf,3):-1:1
    IMt = IMf(:,:,t);
    nanst = isnan(IMt);
    IMt(nanst) = mem(nanst);
    IMf(:,:,t) = gamma*mem + (1-gamma)*IMt;
    mem = IMf(:,:,t);
end
IMf(nans) = nan;

%Difference of Gaussians, normalized over observed pixels: near unscanned
%(NaN) pixels the kernels renormalize instead of averaging in zeros, which
%would attenuate and displace edge pixels by an amount that changes with the
%frame-to-frame NaN pattern. Averaging fewer pixels leaves more noise, so the
%result is divided by its expected noise SD relative to fully observed pixels,
%sqrt(sum(k.^2.*W)/sum(k.*W).^2 / sum(k.^2)) for the center kernel k.
%Identical to plain filtering where all pixels are observed.
chunk = 100; %frames per chunk, to bound memory
for t0 = 1:chunk:nTimePoints
    ix = t0:min(nTimePoints, t0+chunk-1);
    W = single(~nans(:,:,ix));
    X = IMf(:,:,ix); X(nans(:,:,ix)) = 0;
    Wc = imgaussfilt(W, [sigma sigma]);
    G1 = imgaussfilt(X, [sigma sigma]) ./ Wc;
    G1(nans(:,:,ix)) = 0;
    G2 = imgaussfilt(G1, 5*[sigma sigma]) ./ imgaussfilt(W, 5*[sigma sigma]);
    noiseSD = sqrt(imgaussfilt(W, [sigma sigma]/sqrt(2)) ./ Wc.^2); %k.^2 is a Gaussian of width sigma/sqrt(2)
    IMf(:,:,ix) = (G1 - G2) ./ noiseSD;
end
clear W X Wc G1 G2 noiseSD
IMf(nans) = nan;

%nonmax suppression- find maxima. Unobserved pixels are -Inf here so they
%cannot veto a neighbouring or adjacent-frame maximum (NaN would).
%Scores are kept separately by the number of observed neighbours (0-8),
%see below.
nPx = prod(sz(1:2));
nbKernel = ones(3); nbKernel(2,2) = 0;
skByNb = zeros([sz(1:2) 9]); obsByNb = zeros([sz(1:2) 9]);
for fr = size(IMf,3)-ceil(1.5*tau):-1:2 %ceil(tau) because the filtering is uncertain in the final frames

    IMfr = IMf(:,:,fr); IMpre = IMf(:,:,fr-1); IMpost = IMf(:,:,fr+1);
    obs = ~isnan(IMfr);
    IMfr(~obs) = -inf; IMpre(isnan(IMpre)) = -inf; IMpost(isnan(IMpost)) = -inf;

    selMax = IMfr==ordfilt2(IMfr,9, ones(3));
    localMax = selMax & IMfr>IMpre & IMfr>=IMpost;
    %maxinds = find(IMfr==ordfilt2(IMfr,9, ones(3)));
    % sel = IMfr(maxinds)>0 & IMpre(maxinds)<=IMfr(maxinds) & IMpost(maxinds)<=IMfr(maxinds);
    % %sel = IMpre(maxinds)<=IMfr(maxinds) & IMpost(maxinds)<=IMfr(maxinds);
    % maxinds = maxinds(sel);
    nb = conv2(double(obs), nbKernel, 'same'); %observed neighbours
    ix = find(obs);
    ixNb = ix + nb(ix)*nPx;
    obsByNb(ixNb) = obsByNb(ixNb) + 1;
    ix = find(localMax);
    ixNb = ix + nb(ix)*nPx;
    skByNb(ixNb) = skByNb(ixNb) + IMfr(ix).^2;
end
%A pixel with unobserved neighbours (at the edge of the scanned region) has
%fewer neighbours to beat, so it is scored as a maximum more often, which
%raises its background level and lowers that of the pixels inside it. Scale
%each neighbour-count class so that its typical (median over pixels) score
%per observation matches that of fully surrounded pixels; sparse synapses
%hardly move these medians, and signal and background are scaled alike.
classLevel = nan(1,9);
for k = 1:9
    o = obsByNb(:,:,k); s = skByNb(:,:,k);
    enough = o >= 200;
    if nnz(enough) >= 30
        classLevel(k) = median(s(enough)./o(enough));
    end
end
gain = classLevel(9)./classLevel;
gain(~isfinite(gain)) = 1;
skIm = sum(skByNb .* reshape(gain, 1, 1, 9), 3);
nObs = sum(obsByNb, 3);
clear skByNb obsByNb
%unnormalized sums, so images can be accumulated across trials weighted by observations
acc.num = skIm;
acc.count = nObs;
acc.reg = 300;
acc.nFrames = numel(2:size(IMf,3)-ceil(1.5*tau));
skIm = skIm./(acc.reg+nObs); %normalize to # observations w regularizer
clear nans

%summary = skewness(IMf(:,:, 1:end-3*ceil(tau)), 1,3); %.*IMgamma; 
summaryEroded = skIm;

summaryEroded = summaryEroded ./ 10^(floor(log10(max(summaryEroded(:))))-1);
summaryEroded(~valid) = nan;
mfSummary = nanmedfilt2(summaryEroded, [5 5]);
summaryEroded = summaryEroded - mfSummary;
%valid = valid & (skIm ~= 0);
summaryEroded(~valid) = nan;

skIm(~valid) = nan;

thetaf = getActImPeaks(summaryEroded,params.peakth,[],params.minPeakDistance);

P.row = thetaf(:,2);
P.col = thetaf(:,3);
P.val = thetaf(:,1);

peaksMask = zeros(size(summaryEroded));
peaksMask(sub2ind(size(peaksMask), round(P.row), round(P.col))) = 1;

P.peakIM = summaryEroded .* peaksMask;

if doPlot
    figure, imagesc(summaryEroded); %hAx1 = gca;
    hold on, scatter(P.col, P.row, 20*P.val, 'r' ); %'margeredgecolor', 'r');
    figure, imagesc(summaryEroded); %hAx2 = gca;
end

end
