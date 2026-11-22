function [skIm, P] = localizeSources_vIM(IM, vIM, params, doPlot)
%inputs:
%IM:        3D recording, X x Y x Time
%aData:     alignment metadata
%IM and vIM are masked and filtered in place: pass them without keeping
%other references to them (as loadAndProcessTrialAsync does), or MATLAB
%copies both. Movie-sized arrays are cleared as soon as they are used.
tau = params.tau_s.*params.alignHz; %time constant in frames
params.tau_frames = tau;
sigma = params.sigma_px; %space constant in pixels
baselineWindow = ceil(params.baselineWindow_Glu_s.*params.alignHz);
denoiseWindow = ceil(params.denoiseWindow_s.*params.alignHz);
nans = isnan(IM);

sz = size(IM);
valid = mean(nans,3)<params.nanThresh; %a pixel must be imaged at least (1-nanThresh) of the time to be included
clear nans
if ~any(valid)
    warning('Recording had no valid pixels; likely too much motion')
    P = [];
    skIm = nan(sz(1:2));
    return
end
if nargin<4
    doPlot = false;
end

if isempty(vIM) %for raster imaging
    vIM = ones(sz, 'like', IM);
end

%initialize filtered image
IMf = IM; clear IM;
IMf = setFrames(IMf, ~valid, nan);
nans = isnan(IMf);
nanFrac = mean(nans,3);
vIM = setFrames(vIM, nans, nan); %1000*mean(vIM(:,:, 1:min(end,400)), 'all', 'omitnan');

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

%variance of the smoothed data, for Vb and Vk below (computed before the
%baseline so that IMs can be cleared as soon as the baseline is made)
firstValidFrames = find(any(~nans, [1 2]),500, 'first');
varIM = var(IMs(:,:,firstValidFrames),0,3,"omitmissing");

%baseline estimate
IMb = smoothdata(IMs, 3, 'movmedian', baselineWindow, 'omitnan');
clear IMs

%estimate Vb and Vk, parameters for estimating variance from baseline brightness
% Vb: the variance of a 'dim' pixel due to electronic and dark noise
% Vk: the slope of the variance-brightness relationship
varIM(nanFrac>0.4) = nan;
Vb = params.VIF*prctile(varIM, 10, 'all');
varPred = mean(IMb(:,:,firstValidFrames),3,'omitmissing').* mean(vIM(:,:,firstValidFrames),3,'omitmissing');
selBright = varPred>prctile(varPred(:), 90);
Vk = prctile((varIM(selBright)-(Vb/params.VIF))./varPred(selBright), 10);

%Highpass filter in time; This must occur before DoG to avoid edge artifacts
IMf = IMf - IMb; 

stdIM = IMb; clear IMb %in the baseline's memory
stdIM = sqrt(max(0,Vk.*stdIM.*vIM)+Vb); %compute standard deviation
clear vIM
%divide by uncertainty to get a Z-score
IMf = IMf./stdIM;
clear stdIM

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

%Difference of Gaussians
IMf = setFrames(IMf, nans, 0);
IMf = imgaussfilt(IMf, [sigma sigma]);
IMf = IMf - imgaussfilt(IMf, 5*[sigma sigma]);
IMf = setFrames(IMf, nans, nan);

%nonmax suppression- find maxima
skIm = zeros(sz(1:2));
for fr = size(IMf,3)-ceil(1.5*tau):-1:2 %ceil(tau) because the filtering is uncertain in the final frames
    
    IMfr = IMf(:,:,fr); IMpre = IMf(:,:,fr-1); IMpost = IMf(:,:,fr+1);
    
    selMax = IMfr==ordfilt2(IMfr,9, ones(3));
    localMax = selMax & IMfr>IMpre & IMfr>=IMpost;
    %maxinds = find(IMfr==ordfilt2(IMfr,9, ones(3)));
    % sel = IMfr(maxinds)>0 & IMpre(maxinds)<=IMfr(maxinds) & IMpost(maxinds)<=IMfr(maxinds);
    % %sel = IMpre(maxinds)<=IMfr(maxinds) & IMpost(maxinds)<=IMfr(maxinds);
    % maxinds = maxinds(sel);
    skIm(localMax) = skIm(localMax) + IMfr(localMax).^2;
end
clear IMf
skIm = skIm./(300+sum(~nans(:,:,2:end-ceil(1.5*tau)),3)); %normalize to # observations w regularizer
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

function X = setFrames(X, mask, val)
%X(mask) = val, one frame at a time. mask is rows x cols (every frame) or
%rows x cols x frames. Indexing the whole movie with a logical mask makes an
%index array of 8 bytes per masked element, up to twice the movie.
for t = 1:size(X,3)
    x = X(:,:,t);
    x(mask(:,:,min(t,end))) = val;
    X(:,:,t) = x;
end
end
