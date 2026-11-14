function E = processAllTrials_Async(dr, fns, fls, els, selPix, sources, discardFrames, alignData, meanAligned, motOutput, roiData, validTrials, params)
%subproblems run on a pool of worker processes (see extractionPool)
poolTimer = tic;
extractionPool(params);
fprintf('[SILo timing] extraction pool %.1f s\n', toc(poolTimer));
numDatasets = numel(fls);
E = cell(numDatasets,1);
minFrames = 100; %extractTrial initializes from the first 100 samples; shorter trials cannot be source-extracted
pending = []; %trial whose extraction is still running
for i = 1:numel(validTrials)
    nLoad = validTrials(i);
    loadTimer = tic;
    CD = loadTrial(dr, fns{nLoad},fls(nLoad),els(nLoad),selPix,discardFrames{nLoad}, alignData{nLoad}, meanAligned(:,:,:,nLoad), motOutput(:,nLoad), roiData, params);
    E{nLoad}.ROIs = CD.ROIs; E{nLoad}.global = CD.global;
    E{nLoad}.discardFrames = CD.discardFrames;
    E{nLoad}.frameLines = CD.frameLines;
    fprintf('[SILo timing] trial %d: load %.1f s\n', nLoad, toc(loadTimer));
    if ~isempty(pending) %we process the previous trial after loading the next, to keep CPU usage up during loading
        waitTimer = tic;
        E{pending} = processResult(resultsFuture, E{pending},params);
        fprintf('[SILo timing] trial %d: waited %.1f s for extraction\n', pending, toc(waitTimer));
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
    waitTimer = tic;
    E{pending} = processResult(resultsFuture,E{pending},params);
    fprintf('[SILo timing] trial %d: waited %.1f s for extraction\n', pending, toc(waitTimer));
    if isempty(E{pending})
        disp(['Error processing trial: ' int2str(pending) '\n Continuing...'])
    end
end

end
% ---------------- Helper Functions ----------------
function E = processResult(resultsFuture, E, params)
try
    [H,S,LS,F0,SNR] = fetchOutputs(resultsFuture);
    discard = E.discardFrames;
    E.footprints = H;
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
end
%free the trial's futures: kept, the client grew by ~3 GB per 780125 trial
%until the pool was deleted
delete(resultsFuture.Predecessors);
delete(resultsFuture);
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

        alignDataSlap2 = alignData.slap2;

        %upsample motion. If registration undid the online motion correction,
        %each DS frame's motion is relative to the DMD coordinates of its own
        %offset; adding that offset (mapped to rows/cols, relative to the
        %reference) gives the continuous tissue motion, which is what can be
        %interpolated. Frames are then formed compensated, per block, around the
        %block's typical offset (see onlineCompensation).
        [comp, motionDSr, motionDSc] = onlineCompensation(alignDataSlap2, alignData, S2data);
        if numel(alignData.DSframes) > 1
            motionC = interp1(alignData.DSframes, motionDSc, frameLines, 'pchip', 'extrap') + motOutput(2);
            motionR = interp1(alignData.DSframes, motionDSr, frameLines, 'pchip', 'extrap') + motOutput(1);
        else %only one downsampled frame; assume constant motion
            motionC = repmat(motionDSc(1), size(frameLines)) + motOutput(2);
            motionR = repmat(motionDSr(1), size(frameLines)) + motOutput(1);
        end
        CD.motionC = motionC; %with compensation: relative to the frame of comp.RefOffset
        CD.motionR = motionR;

        viewC = alignDataSlap2.viewC(1,:);
        viewR = alignDataSlap2.viewR(:,1);

        selPx2D = selPx2D(1:size(alignDataSlap2.viewC,1), 1:size(alignDataSlap2.viewC,2));
        nPx = numel(selPx2D);
        meanIM = meanIM(1:size(alignDataSlap2.viewC,1), 1:size(alignDataSlap2.viewC,2),:);
        roiData = shiftRoiData(roiData, [0 0], size(selPx2D)); %masks are on SILo's consensus canvas; crop/pad them to this trial's grid, which shares its top left corner
        labeled = medfilt2(meanIM(:,:,params.activityChannel), [3 3]);
        labeled = ~isnan(meanIM(:,:,params.activityChannel)) & labeled>3*prctile(labeled(~isnan(labeled)), 25); %labeled pixels
        
        meanPx = reshape(meanIM, numel(labeled), numChannels);
        meanPx = meanPx(:, orderedChannels); %meanIM is in raw channel order; match the reordered frames below
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
            %and interpolated with interpFrames, at the sampled pixels only, but with
            %the channels in analysis order (getImages returns them in file order)
            [compB, dB] = blockCompensation(comp, frameLines(fIxs));
            [Yq, Fq] = slap2FramesAt(S2data, orderedChannels, frameLines(fIxs), ceil(dt), ...
                alignDataSlap2.trimRows, alignDataSlap2.trimCols, rD, cD, motionR(fIxs) - dB(1), motionC(fIxs) - dB(2), compB);

            IMsel(:, :, fIxs) = Yq(qpos(selPx2D(:)),:,:);
            Finvsel(:,fIxs) = Fq(qpos(selPx2D(:)),:);

            %compute global ROI activity
            yLabeled = double(Yq(qpos(labeled(:)),:,:)); nans= isnan(yLabeled); %per-channel mask, same size as M
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
        meanPx = meanPx(:, orderedChannels); %meanIM is in raw channel order; match the reordered frames below
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
            nanMask = isnan(Y(:,:,orderedChannels,:)); %Y2 below is in analysis channel order
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
            yLabeled = double(Y2(labeled(:),:,:)); nans= isnan(yLabeled); %per-channel mask, same size as M
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

function [Yq, Fq] = slap2FramesAt(S2, ch, times, dt, trimRows, trimCols, rD, cD, mR, mC, comp)
%the frames getImages(ch, times, dt, 1, 1) returns, trimmed to trimRows x
%trimCols and interpolated as interpFrames does, at rows rD + mR(f) and columns
%cD + mC(f) of frame f: Yq is numel(rD) x channels (in the order of ch) x frames, Fq (the variance
%multiplier) numel(rD) x frames. Frames are formed only at the raster
%superpixels that these points' neighbours can fall on in this block, which is
%several times faster than forming whole frames. comp (optional; see
%rawFramesAt) forms the frames with the online motion correction undone, in
%which case mR, mC must be motion estimated on such frames.
if nargin < 11, comp = []; end
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
[V, W] = rawFramesAt(S2, ch, times, dt, id(:), comp);   % needed pixels x frames (x channels)
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

function [comp, motionDSr, motionDSc] = onlineCompensation(alignDataSlap2, alignData, S2)
%the online-compensation settings registration used (empty without), and the
%DS motion to interpolate: as saved, or with compensation the continuous
%motion relative to the frame of the reference offset
motionDSr = double(alignData.motionDSr); motionDSc = double(alignData.motionDSc);
comp = [];
if ~isfield(alignDataSlap2, 'onlineCompensation'), return, end
oc = alignDataSlap2.onlineCompensation;
if ~isfield(oc, 'enabled') || ~logical(oc.enabled), return, end
off = getLineOnlineOffsets(S2, YLead = double(oc.yLead(:))');
map = double(oc.mapping); ref = double(oc.refOffset(:))';
d = map * ([double(oc.frameRefX(:))'; double(oc.frameRefY(:))'] - ref');   %[row; col] per DS frame
motionDSr = reshape(motionDSr, 1, []) + d(1, :);
motionDSc = reshape(motionDSc, 1, []) + d(2, :);
comp = struct('Offsets', off, 'RefOffset', ref, 'Mapping', map);
end

function [compB, dB] = blockCompensation(comp, lines)
%frames of a block are formed compensated around the block's median offset,
%so that the motion to apply, motion - dB, stays as small as without
%compensation; empty and [0 0] without compensation
compB = []; dB = [0 0];
if isempty(comp), return, end
x = double(comp.Offsets.x(lines)); y = double(comp.Offsets.yEff(lines));
oB = round([median(x) median(y)]);
dB = (comp.Mapping * (oB - comp.RefOffset)')';
compB = struct('Offsets', comp.Offsets, 'RefOffset', oB, 'Mapping', comp.Mapping);
end

function [V, W] = rawFramesAt(S2, ch, times, dt, id, comp)
%getImages' frames at raster superpixels id only: V is numel(id) x frames x
%channels (in the order of ch), W (getImages' dataCount) numel(id) x frames.
%As getImages, a frame leaves out the lines exactly dtRead from its time and
%does not skip invalid samples. comp, when not empty, is a struct with fields
%Offsets (getLineOnlineOffsets(S2)) and RefOffset ([x y], or empty for the
%median): the lines' samples are then moved by their online offsets first, and
%id are pixels of the frame RefOffset (see slap2FrameSums).
if nargin < 6, comp = []; end
args = {};
if ~isempty(comp)
    args = {'CompensateOnline', true, 'Offsets', comp.Offsets, 'RefOffset', comp.RefOffset};
    if isfield(comp, 'Mapping'), args = [args {'Mapping', comp.Mapping}]; end
end
[Num, Den, pix] = slap2FrameSums(S2, ch, times, dt, 'Boundary', 'exclusive', 'MaskInvalid', false, ...
    'Pixels', id, args{:});
[tf, loc] = ismember(id(:), pix);
V = nan(numel(id), numel(times), numel(ch));
W = zeros(numel(id), numel(times));
V(tf, :, :) = Num(loc(tf), :, :) ./ Den(loc(tf), :);
W(tf, :) = Den(loc(tf), :);
end
