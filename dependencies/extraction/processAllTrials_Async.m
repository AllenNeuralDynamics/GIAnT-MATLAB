function E = processAllTrials_Async(dr, fns, fls, els, selPix, sources, discardFrames, alignData, meanAligned, motOutput, roiData, validTrials, params)
%subproblems run on a pool of worker processes (see extractionPool)
extractionPool(params);
numDatasets = numel(fls);
E = cell(numDatasets,1);
minFrames = 100; %extractTrial initializes from the first 100 samples; shorter trials cannot be source-extracted
pending = []; %trial whose extraction is still running
for i = 1:numel(validTrials)
    nLoad = validTrials(i);
    CD = loadTrial(dr, fns{nLoad},fls(nLoad),els(nLoad),selPix,discardFrames{nLoad}, alignData{nLoad}, meanAligned(:,:,:,nLoad), motOutput(:,nLoad), roiData, params);
    E{nLoad}.ROIs = CD.ROIs; E{nLoad}.global = CD.global;
    E{nLoad}.discardFrames = CD.discardFrames;
    E{nLoad}.frameLines = CD.frameLines;
    if ~isempty(pending) %we process the previous trial after loading the next, to keep CPU usage up during loading
        [E{pending}, B] = processResult(resultsFuture, E{pending},params);
        if isempty(E{pending})
            disp(['Error processing trial: ' int2str(pending) '\n Continuing...'])
        end
        pending = [];
    end
    nFrames = size(CD.Yobs,3);
    if nFrames < minFrames
        warning('SILo:shortTrial', 'Trial %d (%s) has only %d frames (minimum %d); skipping source extraction for this trial.', ...
            nLoad, fns{nLoad}, nFrames, minFrames);
        clear CD;
        continue
    end
    disp(['Processing trial ' int2str(validTrials(i)) ' from ' fns{nLoad}])
    %Y = squeeze(CD.Yobs(:,1,:));
    Y = permute(CD.Yobs, [1 3 2]);
    resultsFuture = extractTrial(Y,CD.Finv, sources, any(selPix,3), params, alignData{nLoad});
    pending = nLoad;
    clear CD;
end

if ~isempty(pending)
    [E{pending}, B] = processResult(resultsFuture,E{pending},params);
    if isempty(E{pending})
        disp(['Error processing trial: ' int2str(pending) '\n Continuing...'])
    end
end

end
% ---------------- Helper Functions ----------------
function plotE(E, Y,B,selPix) %debugging helper only
    %generate activity movie, baseline movie, and residual
    sel2D = any(selPix,3);
    Ht = reshape(E.footprints,numel(sel2D),[]);
    Ht = Ht(sel2D(:),:);

    Amov = max(0,Ht,'omitmissing')*max(0,E.dF.denoised,'omitmissing');
    Rmov = Y - Amov - B;

    render = zeros([size(sel2D) 500]);
    render(repmat(sel2D,1,1,500)) = B(:,501:1000);
end

function [E, B] = processResult(resultsFuture, E, params)
try
    [H,B,S,LS,F0,SNR] = fetchOutputs(resultsFuture);
    discard = E.discardFrames;
    E.footprints = H;
    %E.baseline = single(B); %This is a huge variable and currently unused.
    E.dF.events = S;
    E.dF.events(:,discard) = nan;

    E.dF.denoised(:,:,1) = convn(S(:,:,1),params.k,'same');
    if size(S,3)>1
        E.dF.denoised(:,:,2) = convn(S(:,:,2),params.k2,'same');
    end
    E.dF.denoised(:,discard,:) = nan;

    E.dF.ls = LS;
    E.dF.ls(:,discard,:) = nan;

    E.F0 = F0;
    E.SNR = SNR;
catch ME
    fprintf(2, 'processResult ERROR: %s\n', ME.message);
    fprintf(2, '%s\n', getReport(ME, 'extended', 'hyperlinks', 'off'));
    B = [];
end
end

function CD = loadTrial(dr, fn, startLine, endLine, selPix, discardFrames, alignData, meanIM, motOutput, roiData, params)
disp('Loading high-res data for file:')
disp([dr filesep fn])

numChannels = params.numChannels;
orderedChannels = [params.activityChannel:numChannels, 1:params.activityChannel-1];

if params.isSLAP2
        if params.includeIntegrationROIs
            warning('includeIntegration not implemented, using raster only!')
            spTypeFlag = 1; %use only raster superpixels
        else
            spTypeFlag = 1; %use only raster superpixels
        end

        %load the high time resolution data
        S2data = slap2.Slap2DataFile([dr filesep fn]);
        meta = loadMetadata([dr filesep fn]);
        
        linerateHz = 1/meta.linePeriod_s;
        dt = linerateHz/params.analyzeHz;
        frameLines = ceil(startLine:dt:endLine);
        nFrames= length(frameLines);
        CD.frameLines = frameLines;
        selPx2D = any(selPix,3);

        %upsample motion
        if numel(alignData.DSframes) > 1
            motionC = interp1(alignData.DSframes, alignData.motionDSc, frameLines, 'pchip', 'extrap') + motOutput(2);
            motionR = interp1(alignData.DSframes, alignData.motionDSr, frameLines, 'pchip', 'extrap') + motOutput(1);
        else %only one downsampled frame; assume constant motion
            motionC = repmat(alignData.motionDSc(1), size(frameLines)) + motOutput(2);
            motionR = repmat(alignData.motionDSr(1), size(frameLines)) + motOutput(1);
        end
        CD.motionC = motionC;
        CD.motionR = motionR;

        alignDataSlap2 = alignData.slap2;
        viewC = alignDataSlap2.viewC(1,:);
        viewR = alignDataSlap2.viewR(:,1);

        selPx2D = selPx2D(1:size(alignDataSlap2.viewC,1), 1:size(alignDataSlap2.viewC,2));
        nPx = numel(selPx2D);
        meanIM = meanIM(1:size(alignDataSlap2.viewC,1), 1:size(alignDataSlap2.viewC,2),:);
        roiData = shiftRoiData(roiData, [0 0], size(selPx2D)); %masks are on SILo's consensus canvas; crop/pad them to this trial's grid, which shares its top left corner
        labeled = medfilt2(meanIM(:,:,params.activityChannel), [3 3]);
        labeled = ~isnan(meanIM(:,:,params.activityChannel)) & labeled>3*prctile(labeled(~isnan(labeled)), 25); %labeled pixels
        
        meanPx = reshape(meanIM, numel(labeled), numChannels);
        mLabeled = meanPx(labeled,:);

        %pixels to sample: sources, labeled pixels (global activity) and user
        %ROIs. Pixel (r,c) of the view is read at row viewR(r) + motionR and
        %column viewC(c) + motionC of the trimmed datafile image.
        need = selPx2D | labeled;
        for rix = 1:length(roiData)
            need = need | reshape(logical(roiData{rix}.mask), size(need));
        end
        q = find(need);
        [rq, cq] = ind2sub(size(need), q);
        rD = viewR(rq); cD = viewC(cq);
        qpos = zeros(nPx, 1); qpos(q) = 1:numel(q); %rows of the sampled pixels q

        blocksize = 600; %number of frames to load at a time
        nBlocks = ceil(nFrames./blocksize);
        blockEdges = round(linspace(1, nFrames+1, nBlocks+1));

        sumF = sum(meanPx,1,'omitmissing');
        IMsel = nan(sum(selPx2D(:)),numChannels, nFrames);
        Finvsel = nan(sum(selPx2D(:)),nFrames);
        for bix = nBlocks:-1:1
            fIxs  = blockEdges(bix):(blockEdges(bix+1)-1);
            nFramesInBlock = length(fIxs);
            %as getImages(orderedChannels, frameLines(fIxs), ceil(dt), 1, 1), trimmed
            %and interpolated with interpFrames, at the sampled pixels only
            [Yq, Fq] = slap2FramesAt(S2data, orderedChannels, frameLines(fIxs), ceil(dt), ...
                alignDataSlap2.trimRows, alignDataSlap2.trimCols, rD, cD, motionR(fIxs), motionC(fIxs));

            IMsel(:, :, fIxs) = Yq(qpos(selPx2D(:)),:,:);
            Finvsel(:,fIxs) = Fq(qpos(selPx2D(:)),:);

            %compute global ROI activity
            yLabeled = double(Yq(qpos(labeled(:)),:,:)); nans= isnan(yLabeled(:,1,:));
            M = repmat(mLabeled,1,1,nFramesInBlock); M(nans) = nan;
            CD.global.F(orderedChannels,fIxs) = (sum(yLabeled,1, 'omitmissing')./sum(M,1, 'omitmissing')).*sumF;

            %compute user ROI activity
            for rix = length(roiData):-1:1
                mask = roiData{rix}.mask;
                tmp1 = Yq(qpos(mask(:)),:,:);  %the data over the ROI pixels
                Fpx{rix}(:,:,fIxs) = tmp1;
                tmp2 = repmat(meanPx(mask(:),:), 1,1,numel(fIxs)); %the mean image over the ROI pixels
                nans= isnan(tmp1) | isnan(tmp2);
                tmp1(nans) = 0;
                tmp2(nans) = 0;
                CD.ROIs.F(rix,:,fIxs) = (sum(tmp1,1)./sum(tmp2,1)).*sum(meanPx(mask(:),:),1,'omitmissing'); %dFF over the valid pixels, times the mean
            end
        end

        %perform SVD on user ROIs to denoise
        CD.ROIs.Fsvd = nan(length(roiData), numChannels, nFrames);
        for rix = 1:length(roiData)
            if ~isempty(Fpx{rix}) && ~all(isnan(Fpx{rix}(:)))
                for cix = 1:numel(orderedChannels)
                    Dtmp = double(permute(Fpx{rix}(:,cix,:), [1 3 2])); %pixels x frames; squeeze would transpose a one-pixel ROI
                    [UU,SS,VV,bg] = nansvd(Dtmp,3, 10, params.nanThresh);
                    roiLikeness = (abs(mean(UU,1, 'omitnan'))./sqrt(mean(UU.^2,1, 'omitnan'))).*diag(SS)'; %components nansvd padded are NaN, which max skips
                    [~,selPC] = max(roiLikeness);
                    CD.ROIs.Fsvd(rix,cix,:) = mean(bg+(UU(:,selPC)*SS(selPC,selPC)*VV(:,selPC)'),1, 'omitnan');
                end
            end
        end

    else
        if endsWith(fn, '.h5')
            desc = h5info([dr filesep fn]);
            IM = h5read([dr filesep fn], ['/', desc.Datasets.Name]);
        else
            IM = ScanImageTiffWrapper([dr filesep fn]);
        end
        IM = double(IM);

        %rearrange into correct dimensions
        selPx2D = any(selPix,3);
        IM = reshape(IM, size(IM,1), size(IM,2), numChannels, []);
        baseline = prctile(reshape(median(IM,4,'omitmissing'), [], numChannels), 5,1);
        IM = IM-reshape(baseline, [1 1 numChannels 1]); %subtract baseline
        nPx = size(IM,1)*size(IM,2);
        nFrames = size(IM,4);
        CD.frameLines = 1:nFrames;

        %upsample motion
        viewC = (1:size(IM,2)) + motOutput(2);
        viewR = (1:size(IM,1))' + motOutput(1);

        CD.motionC = alignData.motionC;
        CD.motionR = alignData.motionR;
        
        meanIM = meanIM(1:size(IM,1), 1:size(IM,2),:);
        roiData = shiftRoiData(roiData, [0 0], [size(IM,1) size(IM,2)]); %masks are on SILo's consensus canvas; crop/pad them to this trial's grid
        scored = medfilt2(meanIM(:,:,params.activityChannel), [3 3]);
        scored = (scored-prctile(scored(:),25,'all'))./(prctile(scored(:),99,'all')-prctile(scored(:),25,'all'));
        labeled = scored>0.1;
        meanPx = reshape(meanIM, numel(labeled), numChannels);
        mLabeled = meanPx(labeled,:);

        blocksize = 600; %number of frames to load at a time;
        nBlocks = ceil(nFrames./blocksize);
        blockEdges = round(linspace(1, nFrames+1, nBlocks+1));

        sumF = sum(meanPx,1,'omitmissing');
        IMsel = nan(sum(selPx2D(:)),numChannels, nFrames);
        Finvsel = nan(sum(selPx2D(:)),nFrames);
        Fresh = ones(size(IM,1), size(IM,2), 'single'); %freshness currently unused for non-SLAP2
        for bix = nBlocks:-1:1
            fIxs  = blockEdges(bix):(blockEdges(bix+1)-1);
            nFramesInBlock = length(fIxs);
            
            Y = IM(:,:,:,fIxs); %reduce communication overhead to parallel workers
            nanMask = isnan(Y);
            Y2 = nan(length(viewR),length(viewC),numChannels, nFramesInBlock);
            Finv = nan(length(viewR),length(viewC),nFramesInBlock);
            parfor frIx = 1:nFramesInBlock
                [Y2(:,:,:,frIx), Finv(:,:,frIx)] = interpFrames(Y(:,:,orderedChannels,frIx),viewC, viewR, Fresh);
            end
            Y2(nanMask) = nan;

            Y2 = reshape(Y2, nPx, numChannels, nFramesInBlock);
            Finv= reshape(Finv, nPx, nFramesInBlock);

            IMsel(:, :, fIxs) = Y2(selPx2D,:,:);
            Finvsel(:,fIxs) = Finv(selPx2D,:);

            %compute global ROI activity
            yLabeled = double(Y2(labeled(:),:,:)); nans= isnan(yLabeled(:,1,:));
            M = repmat(mLabeled,1,1,nFramesInBlock); M(nans) = nan;
            CD.global.F(orderedChannels,fIxs) = (sum(yLabeled,1, 'omitmissing')./sum(M,1, 'omitmissing')).*sumF;

            %compute user ROI activity
            for rix = length(roiData):-1:1
                mask = roiData{rix}.mask;
                tmp1 = Y2(mask(:),:,:);  %the data over the ROI pixels
                Fpx{rix}(:,:,fIxs) = tmp1;
                tmp2 = repmat(meanPx(mask(:),:), 1,1,numel(fIxs)); %the mean image over the ROI pixels
                nans= isnan(tmp1) | isnan(tmp2);
                tmp1(nans) = 0;
                tmp2(nans) = 0;
                CD.ROIs.F(rix,:,fIxs) = (sum(tmp1,1)./sum(tmp2,1)).*sum(meanPx(mask(:),:),1,'omitmissing'); %dFF over the valid pixels, times the mean
            end
        end

        %perform SVD on user ROIs to denoise
        CD.ROIs.Fsvd = nan(length(roiData), numChannels, nFrames);
        for rix = 1:length(roiData)
            if ~isempty(Fpx{rix}) && ~all(isnan(Fpx{rix}(:)))
                for cix = 1:numel(orderedChannels)
                    Dtmp = double(permute(Fpx{rix}(:,cix,:), [1 3 2])); %pixels x frames; squeeze would transpose a one-pixel ROI
                    [UU,SS,VV,bg] = nansvd(Dtmp,3, 10, params.nanThresh);
                    roiLikeness = (abs(mean(UU,1, 'omitnan'))./sqrt(mean(UU.^2,1, 'omitnan'))).*diag(SS)'; %components nansvd padded are NaN, which max skips
                    [~,selPC] = max(roiLikeness);
                    CD.ROIs.Fsvd(rix,cix,:) = mean(bg+(UU(:,selPC)*SS(selPC,selPC)*VV(:,selPC)'),1, 'omitnan');
                end
            end
        end
end

if numel(discardFrames) > 1
    discard = interp1(1:numel(discardFrames), double(discardFrames(:)), linspace(1, numel(discardFrames), size(IMsel,3)))>0; %upsample the discard frames
else
    discard = repmat(any(discardFrames(:)), 1, size(IMsel,3));
end
IMsel(:,:,discard) = nan;     %throw away movement frames as above
IMsel(IMsel<-2*std(IMsel,0,3, 'omitmissing')) = nan; %throw away any aberrantly negative data
CD.global.F(:,discard) = nan;
CD.ROIs.F(:, :, discard) = nan;
CD.ROIs.Fsvd(:,:,discard) = nan;

% Remove nans from IMsel and Finvsel
nans = isnan(IMsel);
nanFill = IMsel;
nanFill(all(nans,3)) = 0;
while any(isnan(nanFill), 'all')
    nanFill = smoothdata(nanFill,3,"movmean",params.baselineWindow_samps, 'omitmissing');
end
IMsel(nans) = nanFill(nans); 
Finvsel(squeeze(nans(:,1,:))) = 1000*mean(Finvsel,'all', 'omitmissing');

CD.Yobs = IMsel;
CD.Finv = Finvsel;
CD.discardFrames = discard;
end

function [Yq, Fq] = slap2FramesAt(S2, ch, times, dt, trimRows, trimCols, rD, cD, mR, mC)
%the frames getImages(ch, times, dt, 1, 1) returns, trimmed to trimRows x
%trimCols and interpolated as interpFrames does, at rows rD + mR(f) and columns
%cD + mC(f) of frame f: Yq is numel(rD) x channels x frames, Fq (the variance
%multiplier) numel(rD) x frames. Frames are formed only at the raster
%superpixels that these points' neighbours can fall on in this block, and only
%the lines that sample those are read, which is several times faster than
%forming whole frames.
hdr = S2.hMultiDataFiles.header;
NR = double(hdr.dmdPixelsPerRow); NC = double(hdr.dmdPixelsPerColumn);
nF = numel(times); nCh = numel(ch); nq = numel(rD); n1 = numel(trimRows); n2 = numel(trimCols);
%pixels within the block's motion range of each point
mark = false(n1, n2);
r0 = floor(rD(:) + min(mR)); c0 = floor(cD(:) + min(mC));
for dr = 0:(ceil(max(mR) - min(mR)) + 1)
    for dc = 0:(ceil(max(mC) - min(mC)) + 1)
        r = r0 + dr; c = c0 + dc; ok = r >= 1 & r <= n1 & c >= 1 & c <= n2;
        mark(r(ok) + n1 * (c(ok) - 1)) = true;
    end
end
need = find(mark); [rt, ct] = ind2sub([n1 n2], need);
id = trimCols(ct(:)) + NR * (trimRows(rt(:)) - 1);   % getImages' raster image is the transpose of superpixel order
[V, W] = rawFramesAt(S2, ch, times, dt, id(:), NR * NC);   % needed pixels x frames (x channels)
nN = numel(need);
V(nN+1, :, :) = nan; W(nN+1, :) = nan;   % the row for neighbours outside the image
tpos = (nN + 1) * ones(n1 * n2 + 1, 1); tpos(need + 1) = 1:nN;   % trimmed-image index + 1 -> row
%the four neighbours of each point in each frame, and their bilinear weights (as interpFrames)
R = rD(:) + mR(:).'; C = cD(:) + mC(:).';
rf = floor(R); cf = floor(C); fr = R - rf; fc = C - cf;
rc = ceil(R); ccl = ceil(C);
fOff = 1:nF;   % frame index, along columns
num = zeros(nq, nF, nCh); den = zeros(nq, nF);
rs = {rf, rf, rc, rc}; cs = {cf, ccl, cf, ccl};
for k = 1:4
    switch k
        case 1, w = (1 - fc) .* (1 - fr);
        case 2, w = fc .* (1 - fr);
        case 3, w = (1 - fc) .* fr;
        case 4, w = fc .* fr;
    end
    r = rs{k}; c = cs{k};
    ok = r >= 1 & r <= n1 & c >= 1 & c <= n2;
    j = tpos((r + n1 * (c - 1)) .* ok + 1);    % invalid -> the NaN column
    ix = j + (nN + 1) * (fOff - 1);   % neighbouring points read neighbouring rows
    f = W(ix);
    den = den + w .* f;
    wf = w .* f;
    for cI = 1:nCh
        num(:, :, cI) = num(:, :, cI) + wf .* V(ix + (nN + 1) * nF * (cI - 1));
    end
end
Yq = permute(num ./ den, [1 3 2]);
Fq = 1 ./ den; Fq(isinf(Fq)) = nan;
end

function [V, W] = rawFramesAt(S2, ch, times, dt, id, nRaster)
%getImages' frames at raster superpixels id only: V is numel(id) x frames x
%channels, W (getImages' dataCount) numel(id) x frames
hdr = S2.hMultiDataFiles.header; P = double(hdr.linesPerCycle);
dtRead = max(3 * dt, P); total = S2.hMultiDataFiles.totalNumLines;
inWin = false(1, max(times) + dtRead);
for t = times(:).', inWin(max(1, t - dtRead):min(total, t + dtRead)) = true; end
lineTimes = find(inWin);
li = mod(lineTimes - 1, P) + 1; ci = floor((lineTimes - 1) / P) + 1;
z = S2.hMultiDataFiles.lineFastZIdxs(li) == 1;
lineTimes = lineTimes(z); li = li(z); ci = ci(z);
map = zeros(nRaster, 1); map(id) = 1:numel(id);
nF = numel(times); nN = numel(id); nCh = numel(ch);
Ls = unique(li); sel = cell(1, P); col = cell(1, P); use = false(size(Ls));
for j = 1:numel(Ls)
    ids = double(S2.hMultiDataFiles.lineSuperPixelIDs{Ls(j)}(:)) + 1;
    m = zeros(size(ids)); ok = ids <= nRaster; m(ok) = map(ids(ok));
    sel{Ls(j)} = find(m > 0); col{Ls(j)} = m(m > 0); use(j) = any(m > 0);
end
keep = ismember(li, Ls(use));   % only lines that sample a needed superpixel are read
lineTimes = lineTimes(keep); li = li(keep); ci = ci(keep);
ld = S2.hMultiDataFiles.getLineData(li, ci, ch);
Lu = Ls(use); parts = cell(numel(Lu), nCh); wparts = cell(numel(Lu), 1); cols = cell(numel(Lu), 1);
for j = 1:numel(Lu)
    L = Lu(j); k = find(li == L);
    D = cat(3, ld{k}); D = double(D(sel{L}, :, :));   % pixels x channels x cycles
    d = times(:) - lineTimes(k);
    [fi, ki] = find(abs(d) < dtRead);
    Wk = sparse(fi, ki, exp(-abs(d(abs(d) < dtRead)) / dt), nF, numel(k));   % frames x cycles
    for c = 1:nCh, parts{j, c} = Wk * reshape(D(:, c, :), numel(sel{L}), []).'; end
    wparts{j} = full(sum(Wk, 2)) * ones(1, numel(sel{L}));
    cols{j} = col{L};
end
cols = vertcat(cols{:}); S = sparse(1:numel(cols), cols, 1, numel(cols), nN);   % sums pixels sampled on several lines
V = zeros(nN, nF, nCh);
for c = 1:nCh, V(:, :, c) = full(S.' * [parts{:, c}].'); end
W = full(S.' * [wparts{:}].');
V = V ./ W;
end
