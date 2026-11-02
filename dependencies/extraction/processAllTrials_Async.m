function E = processAllTrials_Async(dr, fns, fls, els, selPix, sources, discardFrames, alignData, meanAligned, warp, roiData, validTrials, params)
%warp: the across-trial alignment from SILo, mapping the common (across-trial
%aligned) frame to each trial's registered image (see trialWarpCoords)
curPool = gcp('nocreate');
if isempty(curPool) || ~strcmpi(class(curPool), 'parallel.ThreadPool') %  ~strcmpi(class(curPool), 'parallel.ProcessPool') %
    delete(curPool);
    parpool('Threads'); %use all available Threads
end
numDatasets = numel(fls);
E = cell(numDatasets,1); 
for i = 1:numel(validTrials)
    nLoad = validTrials(i);
    CD = loadTrial(dr, fns{nLoad},fls(nLoad),els(nLoad),selPix,discardFrames{nLoad}, alignData{nLoad}, meanAligned(:,:,:,nLoad), warp, nLoad, roiData, params);
    E{nLoad}.ROIs = CD.ROIs; E{nLoad}.global = CD.global;
    E{nLoad}.discardFrames = CD.discardFrames;
    E{nLoad}.frameLines = CD.frameLines;
    if i>1 %we process the previous trial after loading the next, to keep CPU usage up during loading
        [E{validTrials(i-1)}, B] = processResult(resultsFuture, E{validTrials(i-1)},params);
        if isempty(E{nLoad})
            disp(['Error processing trial: ' int2str(nLoad) '\n Continuing...'])
        end
    end
    disp(['Processing trial ' int2str(validTrials(i)) ' from ' fns{nLoad}])
    %Y = squeeze(CD.Yobs(:,1,:));
    Y = permute(CD.Yobs, [1 3 2]);
    prm = params;
    prm.noiseACF = CD.noiseACF;
    if any(CD.noiseACF)
        disp(['Noise correlation between frames at lags 1-' int2str(numel(CD.noiseACF)) ': ' mat2str(CD.noiseACF, 2)])
    end
    resultsFuture = extractTrial(Y,CD.Finv, sources, any(selPix,3), prm, alignData{nLoad});
    clear CD;
end

[E{nLoad}, B] = processResult(resultsFuture,E{nLoad},params);
if isempty(E{nLoad})
    disp(['Error processing trial: ' int2str(nLoad) '\n Continuing...'])
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

    %per-pixel background averaged over kept frames (selected pixels x
    %analyzed channels); a proxy for resting brightness, e.g. spine size.
    %Frames near trial edges are left out, where the baseline fit is biased
    %(GitHub issue #17); trials too short to trim are averaged whole.
    keep = ~discard;
    nTrim = round(params.F0imTrim_s * params.analyzeHz);
    trimmed = keep;
    trimmed([1:min(nTrim,end), max(1,end-nTrim+1):end]) = false;
    if any(trimmed)
        keep = trimmed;
    end
    E.F0im = single(reshape(mean(B(:,keep,:),2,'omitmissing'), size(B,1), size(B,3)));
catch ME
    fprintf(2, 'processResult ERROR: %s\n', ME.message);
    fprintf(2, '%s\n', getReport(ME, 'extended', 'hyperlinks', 'off'));
    B = [];
end
end

function CD = loadTrial(dr, fn, startLine, endLine, selPix, discardFrames, alignData, meanIM, warp, trialIx, roiData, params)
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

        %upsample motion (within-trial, rigid)
        motionC = interp1(alignData.DSframes, alignData.motionDSc, frameLines, 'pchip', 'extrap');
        motionR = interp1(alignData.DSframes, alignData.motionDSr, frameLines, 'pchip', 'extrap');

        alignDataSlap2 = alignData.slap2;
        viewC = alignDataSlap2.viewC(1,:);
        viewR = alignDataSlap2.viewR(:,1);
        assert(all(diff(viewC)==1) && all(diff(viewR)==1), 'loadTrial:view', 'expected unit-step viewC/viewR');

        %record motion including the across-trial offset at the image centre
        x0 = [(numel(viewR)+1)/2; (numel(viewC)+1)/2];
        [rC, cC] = trialWarpCoords(warp, trialIx, x0(1), x0(2));
        dCentre = [rC; cC] - x0;
        CD.motionC = motionC + dCentre(2);
        CD.motionR = motionR + dCentre(1);

        selPx2D = selPx2D(1:size(alignDataSlap2.viewC,1), 1:size(alignDataSlap2.viewC,2));
        nPx = numel(selPx2D);
        meanIM = meanIM(1:size(alignDataSlap2.viewC,1), 1:size(alignDataSlap2.viewC,2),:);
        labeled = medfilt2(meanIM(:,:,params.activityChannel), [3 3]);
        labeled = ~isnan(meanIM(:,:,params.activityChannel)) & labeled>3*prctile(labeled(~isnan(labeled)), 25); %labeled pixels
        
        meanPx = reshape(meanIM, numel(labeled), numChannels);
        mLabeled = meanPx(labeled,:);

        %pixels to sample: sources, labeled pixels (global activity) and user
        %ROIs. A common-frame pixel sits at trialWarpCoords(warp, trialIx, row, col)
        %in this trial's registered image, which is offset from the datafile by viewR/viewC.
        need = selPx2D | labeled;
        for rix = 1:length(roiData)
            if numel(roiData{rix}.mask) == nPx
                need = need | reshape(logical(roiData{rix}.mask), size(need));
            else
                need(:) = true;
            end
        end
        q = find(need);
        [rq, cq] = ind2sub(size(need), q);
        [rT, cT] = trialWarpCoords(warp, trialIx, rq, cq);
        rD = viewR(1) - 1 + rT;
        cD = viewC(1) - 1 + cT;

        %correlation of the noise between analysis frames at the source
        %pixels, which share raw samples (see frameNoiseACF)
        hdr = S2data.hMultiDataFiles.header;
        sel = selPx2D(q);
        noiseMask = []; %all imaged pixels
        if any(sel)
            rM = alignDataSlap2.trimRows(min(max(round(rD(sel)), 1), numel(alignDataSlap2.trimRows)));
            cM = alignDataSlap2.trimCols(min(max(round(cD(sel)), 1), numel(alignDataSlap2.trimCols)));
            noiseMask = false(double(hdr.dmdPixelsPerColumn), double(hdr.dmdPixelsPerRow)); %getImages' raster image size
            noiseMask(sub2ind(size(noiseMask), rM(:), cM(:))) = true;
        end
        nLags = ceil(2 * max(3 * ceil(dt), double(hdr.linesPerCycle)) / dt); %frames sharing raw samples
        CD.noiseACF = frameNoiseACF(S2data, frameLines, ceil(dt), 1, nLags, noiseMask);

        blocksize = 600; %number of frames to load at a time; 100 frames ~= 1GB RAM usage
        nBlocks = ceil(nFrames./blocksize);
        blockEdges = round(linspace(1, nFrames+1, nBlocks+1));

        sumF = sum(meanPx,1,'omitmissing');
        IMsel = nan(sum(selPx2D(:)),numChannels, nFrames);
        Finvsel = nan(sum(selPx2D(:)),nFrames);
        for bix = nBlocks:-1:1
            fIxs  = blockEdges(bix):(blockEdges(bix+1)-1);
            [Y, Fresh] = S2data.getImages(orderedChannels,frameLines(fIxs),ceil(dt),1,1);
            Y = Y(alignDataSlap2.trimRows, alignDataSlap2.trimCols,:,:);
            Fresh = Fresh(alignDataSlap2.trimRows, alignDataSlap2.trimCols,:);
            nFramesInBlock = size(Y,4);

            Yq = nan(numel(q), numChannels, nFramesInBlock);
            Fq = nan(numel(q), nFramesInBlock);
            parfor frIx = 1:nFramesInBlock
                [Yq(:,:,frIx), Fq(:,frIx)] = interpFramesAt(Y(:,:,:,frIx), rD+motionR(fIxs(frIx)), cD+motionC(fIxs(frIx)), Fresh(:,:,frIx));
            end
            Y2 = nan(nPx, numChannels, nFramesInBlock);
            Finv = nan(nPx, nFramesInBlock);
            Y2(q,:,:) = Yq;
            Finv(q,:) = Fq;

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
                    Dtmp = squeeze(double(Fpx{rix}(:,cix,:)));
                    [UU,SS,VV,bg] = nansvd(Dtmp,3, 10, params.nanThresh);
                    roiLikeness = (abs(mean(UU,1, 'omitnan'))./sqrt(mean(UU.^2,1, 'omitnan')))*SS;
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
        CD.noiseACF = 0; %frames are separate exposures

        %upsample motion
        %position of every common-frame pixel in this trial's registered movie
        [cg, rg] = meshgrid(1:size(IM,2), 1:size(IM,1));
        [rq, cq] = trialWarpCoords(warp, trialIx, rg(:), cg(:));
        clear cg rg

        CD.motionC = alignData.motionC;
        CD.motionR = alignData.motionR;
        
        meanIM = meanIM(1:size(IM,1), 1:size(IM,2),:);
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
            Y2 = nan(size(IM,1),size(IM,2),numChannels, nFramesInBlock);
            Finv = nan(size(IM,1),size(IM,2),nFramesInBlock);
            parfor frIx = 1:nFramesInBlock
                [Yf, Ff] = interpFramesAt(Y(:,:,orderedChannels,frIx), rq, cq, Fresh);
                Y2(:,:,:,frIx) = reshape(Yf, size(Fresh,1), size(Fresh,2), []);
                Finv(:,:,frIx) = reshape(Ff, size(Fresh,1), size(Fresh,2));
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
                    Dtmp = squeeze(double(Fpx{rix}(:,cix,:)));
                    [UU,SS,VV,bg] = nansvd(Dtmp,3, 10, params.nanThresh);
                    roiLikeness = (abs(mean(UU,1, 'omitnan'))./sqrt(mean(UU.^2,1, 'omitnan')))*SS;
                    [~,selPC] = max(roiLikeness);
                    CD.ROIs.Fsvd(rix,cix,:) = mean(bg+(UU(:,selPC)*SS(selPC,selPC)*VV(:,selPC)'),1, 'omitnan');
                end
            end
        end
end

discard = interp1(1:numel(discardFrames), double(discardFrames(:)), linspace(1, numel(discardFrames), size(IMsel,3)))>0; %upsample the discard frames
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