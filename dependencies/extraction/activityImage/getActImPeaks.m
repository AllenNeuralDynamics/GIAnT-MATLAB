function thetaf = getActImPeaks(actIM, peakth, exclusionMask, minPeakDistance)
% Defaults
if nargin < 3 || isempty(exclusionMask)
    exclusionMask = false(size(actIM));
end
if nargin < 4 || isempty(minPeakDistance)
    minPeakDistance = 1;
end
minPeakDistance = max(1, minPeakDistance);
peakExclusionSE = ones(2*minPeakDistance - 1);

peakFunc = @gaussianPeaksIntegrated;
ampScale = 1 ./ 0.75;

mu_bg = median(actIM,'all','omitmissing');
sigma_bg = mad(actIM(~isnan(actIM)),1,'all') ./ 0.6741891400433162;
peak_thresh = mu_bg + peakth * sigma_bg;

opts = optimset('MaxFunEvals',5000,'Display','off');

explored = actIM .* ~exclusionMask;
explored(isnan(explored)) = -inf; %ordfilt2 returns NaN for any window touching a NaN, which would rule out peaks next to unscanned pixels
pTmp = ordfilt2(explored, 8, ones(3)) > peak_thresh & ...
       explored == ordfilt2(explored, 9, ones(3));
thetaf = zeros(0,4);
pLocs = zeros(0,2);

if sum(pTmp(:))
    nPeaks = sum(pTmp(:));

    [pY, pX] = ind2sub(size(pTmp),find(pTmp));
    amp = actIM(pTmp) .* ampScale;
    widths = 0.35 * ones(nPeaks,1);
    
    actSelPix = imdilate(pTmp, ones(9)) & ~isnan(actIM);
    [actSelY, actSelX] = ind2sub(size(pTmp),find(actSelPix));
    
    theta0 = [thetaf; [amp, pY, pX, widths]];
    pLocs = [pY, pX];

    ub = [Inf*ones(size(theta0,1),1),min(size(actIM,1)+0.5,pLocs(:,1)+1.5),min(size(actIM,2)+0.5,pLocs(:,2)+1.5),5*ones(size(theta0,1),1)];
    lb = [zeros(size(theta0,1),1),max(0.5,pLocs(:,1)-1.5),max(0.5,pLocs(:,2)-1.5),zeros(size(theta0,1),1)];

    thetaf = lsqcurvefit(peakFunc,theta0,[actSelY,actSelX],actIM(actSelPix)-mu_bg,lb,ub,opts);

    pIM = false(size(actIM));
    iy = min(max(round(thetaf(:,2)), 1), size(actIM,1));
    ix = min(max(round(thetaf(:,3)), 1), size(actIM,2));
    pIM(sub2ind(size(actIM), iy, ix)) = true;
    bufferMask = imdilate(pIM, peakExclusionSE);

    fitIM = zeros(size(actIM));
    fitIM(actSelPix) = peakFunc(thetaf,[actSelY,actSelX]);
    resIM = (actIM - fitIM - mu_bg) ./ sigma_bg;

    fitSupport = (fitIM > 1e-3);
    rejectMask = false(size(actIM));
    explored = resIM .* ~bufferMask .* ~exclusionMask .* ~rejectMask .* fitSupport .* actSelPix;
    pTmp = zeros(size(explored));
    if any(fitSupport(:)) && max(explored(:)) > peakth
        pTmp = explored == max(explored(:));
    end

    while sum(pTmp(:))
        idxNew = find(pTmp,1,'first');
        [pY, pX] = ind2sub(size(pTmp),idxNew);
        amp = actIM(idxNew) .* ampScale;

        nBefore = size(thetaf,1);
        thetaf = [thetaf; [amp, pY, pX, 0.35]];
        pLocs = [pLocs; [pY, pX]];
        newIdx = nBefore + 1;

        CC = bwconncomp(actSelPix);
        foundCC = false;
        for i = 1:CC.NumObjects
            if ismember(idxNew, CC.PixelIdxList{i})
                actSelPix = false(size(actSelPix));
                actSelPix(CC.PixelIdxList{i}) = true;
                [actSelY, actSelX] = ind2sub(size(pTmp),find(actSelPix));
                
                % Find which theta indices correspond to peaks in this connected component
                iy = min(max(round(thetaf(:,2)), 1), size(actIM,1));
                ix = min(max(round(thetaf(:,3)), 1), size(actIM,2));
                peakIndices = sub2ind(size(actIM), iy, ix);
                thetaIdxsToFit = actSelPix(peakIndices);

                foundCC = true;
                break;
            end
        end
        if ~foundCC
            error('getActImPeaks:PeakNotInConnectedComponent', ...
                'Peak linear index %d (row %d, col %d) is not in any connected component of actSelPix (%d objects).', ...
                idxNew, pY, pX, CC.NumObjects);
        end

        ub = [Inf*ones(sum(thetaIdxsToFit),1),min(size(actIM,1)+0.5,pLocs(thetaIdxsToFit,1)+1.5),min(size(actIM,2)+0.5,pLocs(thetaIdxsToFit,2)+1.5),5*ones(sum(thetaIdxsToFit),1)];
        lb = [zeros(sum(thetaIdxsToFit),1),max(0.5,pLocs(thetaIdxsToFit,1)-1.5),max(0.5,pLocs(thetaIdxsToFit,2)-1.5),zeros(sum(thetaIdxsToFit),1)];
        thetaf(thetaIdxsToFit,:) = lsqcurvefit(peakFunc,thetaf(thetaIdxsToFit,:),[actSelY,actSelX],actIM(actSelPix)-mu_bg,lb,ub,opts);

        mu_y_new = thetaf(newIdx,2);
        mu_x_new = thetaf(newIdx,3);
        if abs(mu_y_new - round(mu_y_new)) < 1e-3 && abs(mu_x_new - round(mu_x_new)) < 1e-3
            fitMask = thetaIdxsToFit;
            fitMask(newIdx) = false;
            if any(fitMask)
                ub = [Inf*ones(sum(fitMask),1),min(size(actIM,1)+0.5,pLocs(fitMask,1)+1.5),min(size(actIM,2)+0.5,pLocs(fitMask,2)+1.5),5*ones(sum(fitMask),1)];
                lb = [zeros(sum(fitMask),1),max(0.5,pLocs(fitMask,1)-1.5),max(0.5,pLocs(fitMask,2)-1.5),zeros(sum(fitMask),1)];
                thetaf(fitMask,:) = lsqcurvefit(peakFunc,thetaf(fitMask,:),[actSelY,actSelX],actIM(actSelPix)-mu_bg,lb,ub,opts);
            end
            thetaf(newIdx,:) = [];
            pLocs(newIdx,:) = [];
            rejectMask(idxNew) = true;
        end
        
        pIM = false(size(actIM));
        iy = min(max(round(thetaf(:,2)), 1), size(actIM,1));
        ix = min(max(round(thetaf(:,3)), 1), size(actIM,2));
        pIM(sub2ind(size(actIM), iy, ix)) = true;
        bufferMask = imdilate(pIM, peakExclusionSE);

        fitIM(actSelPix) = peakFunc(thetaf,[actSelY,actSelX]);
        resIM = (actIM - fitIM - mu_bg) ./ sigma_bg;
    
        actSelPix = imdilate(pIM, ones(9)) & ~isnan(actIM);

        fitSupport = (fitIM > 1e-3);
        explored = resIM .* ~bufferMask .* ~exclusionMask .* ~rejectMask .* fitSupport .* actSelPix;
        pTmp = zeros(size(explored));
        if any(fitSupport(:)) && max(explored(:)) > peakth
            pTmp = explored == max(explored(:));
        end
    end

    adj_thresh = peak_thresh ./ (pi/2.*thetaf(:,4).^2.*erf(1./(sqrt(2).*thetaf(:,4))).^2);
    small_peaks = thetaf(:,1) < adj_thresh;
    
    thetaf(small_peaks,:) = [];
end

end

function val = gaussianPeaksIntegrated(theta, yxdata)
% Same interface as original:
% theta  : N×4  [A, mu_y, mu_x, sigma]
% yxdata : M×2  pixel centers [y, x]
% val    : M×1  integrated Gaussian over each pixel

y = yxdata(:,1);
x = yxdata(:,2);

% --- infer pixel size from grid ---
xu = unique(x);
yu = unique(y);

dx = median(diff(xu));
dy = median(diff(yu));

% build edges
xEdges = [xu(1)-dx/2; (xu(1:end-1)+xu(2:end))/2; xu(end)+dx/2];
yEdges = [yu(1)-dy/2; (yu(1:end-1)+yu(2:end))/2; yu(end)+dy/2];

W = numel(xu);
H = numel(yu);

% reshape index map
[~, xIdx] = ismember(x, xu);
[~, yIdx] = ismember(y, yu);

% --- parameters ---
A  = theta(:,1).';   % 1×N
my = theta(:,2).';
mx = theta(:,3).';
s  = max(theta(:,4).', eps);

c   = sqrt(pi/2);
rt2 = sqrt(2);

% --- integrate in x (W×N) ---
xL = xEdges(1:end-1);
xR = xEdges(2:end);
Ix = c .* s .* ( ...
    erf((xR - mx)./(rt2*s)) - ...
    erf((xL - mx)./(rt2*s)) );

% --- integrate in y (H×N) ---
yB = yEdges(1:end-1);
yT = yEdges(2:end);
Iy = c .* s .* ( ...
    erf((yT - my)./(rt2*s)) - ...
    erf((yB - my)./(rt2*s)) );

% --- combine to full image ---
img = (Iy .* A) * Ix.';   % H×W

% --- return in original point ordering ---
val = img(sub2ind([H W], yIdx, xIdx));
end
