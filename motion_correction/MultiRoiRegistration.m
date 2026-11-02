function params = MultiRoiRegistration(fullPathToTrialTable, paramsIn)

if ~nargin
    [fn, trialtabledr] = uigetfile('*.h5', 'Select a trial_table file', '*trial_table*.h5' );
else
    [trialtabledr, fn, ext] = fileparts(fullPathToTrialTable); fn = [fn ext];
end

%PARAMETER SETTING
if nargin>1
    params = setParams('MultiRoiRegistration', paramsIn);
else
    params = setParams('MultiRoiRegistration');
end

params.startTime = char(datetime('now','TimeZone','local','Format','yyyy-MM-dd''T''HH:mm:ss.SSSZZZZZ'));

% Weighted-xcorr acceleration is optional and must never make registration
% dependent on a compiler or MEX binary. Production runs never compile here.
if params.useFastWeightedXcorr && params.useMexWeightedXcorr
    if exist('xcorr2_nans_weighted_mex','file') == 3
        fprintf('Weighted xcorr: optional MEX backend found; each worker will validate it before first use.\n');
    else
        fprintf(['Weighted xcorr: MEX backend not found/loadable; using current MATLAB fast fallback. ' ...
            'Registration will continue normally.\n']);
    end
end
if params.useAdaptiveWeightedXcorr
    warning('MultiRoiRegistration:AdaptiveXcorrEnabled', ...
        ['Adaptive weighted xcorr is enabled. Initial-template correlation remains exhaustive; ' ...
         'main-loop frames use radius %d with full +/- %d fallback/audits.'], ...
        params.adaptiveXcorrRadius,params.clipShift);
end

%load the trial Table, which sets correspondences between the two DMDs
trialTable = loadStructFromH5([trialtabledr filesep fn]);

mocosavedr = fullfile(trialTable.savedr,'motion_correction');
if ~exist(mocosavedr,'dir')
    mkdir(mocosavedr)
end

%set up parallelization
% Use exactly the requested number of process workers (subject to the
% number of trials, logical cores, and the local Processes profile limit).
% The previous implementation could silently use MORE workers than
% params.nWorkers when a larger pool already existed, causing unexpected RAM use.
nDMDs = size(trialTable.filename,1);
nTrials = size(trialTable.true_trial_ix,2);

core_info = evalc('feature(''numcores'');');
core_match = regexp(core_info,'assigned: \d+ logical cores','match','once');
if isempty(core_match)
    numLogicalCores = feature('numcores');
else
    numLogicalCores = sscanf(core_match, 'assigned: %d logical cores');
end

profileMaxWorkers = inf;
try
    localCluster = parcluster('Processes');
    profileMaxWorkers = localCluster.NumWorkers;
catch
    % If the profile cannot be queried, parpool will provide a useful error.
end

requestedWorkers = params.nWorkers;
nWorkers = max(1, min([requestedWorkers, nTrials, numLogicalCores, profileMaxWorkers]));
params.nWorkersRequested = requestedWorkers;
params.nWorkers = nWorkers; % record the number actually used

if nWorkers < requestedWorkers
    warning('MultiRoiRegistration:WorkerLimit', ...
        'Requested %d workers; using %d (limited by trials, logical cores, or Processes profile).', ...
        requestedWorkers, nWorkers);
end

p = gcp('nocreate');
if nWorkers > 1
    if isempty(p) || ~isa(p, 'parallel.ProcessPool') || p.NumWorkers ~= nWorkers
        delete(gcp('nocreate'));
        p = parpool('processes', nWorkers);
    end
    % Workers must resolve GIAnT functions from THIS checkout and must not
    % keep compiled copies from before the current code was pulled. Neither
    % is guaranteed by a reused (or even a fresh) process pool.
    refreshPoolWorkerCode(p, fileparts(fileparts(mfilename('fullpath'))));
else
    % A previously opened process pool still consumes RAM even though the
    % sequential branch below would not use it.
    delete(gcp('nocreate'));
end

trialTable.motion_correction.registration_failed = false(nDMDs, nTrials);

for DMD_ix = 1:nDMDs
    fnRegDS = cell(1,nTrials);
    fnAdata = cell(1,nTrials);
    firstLine = nan(1,nTrials);
    regFail = false(1, nTrials);

    % Broadcast only the fields alignment actually needs. In particular,
    % avoid copying the often very large reference stacks to every worker
    % unless refStackTemplate is enabled.
    workerTable = makeWorkerTable(trialTable, params, DMD_ix);

    if nWorkers>1
        parfor f_ix = 1:nTrials
            [fnRegDS{f_ix}, fnAdata{f_ix}, firstLine(f_ix), regFail(f_ix)]= alignAsync(workerTable, params, f_ix, DMD_ix);
        end
    else
        for f_ix = 1:nTrials
            [fnRegDS{f_ix}, fnAdata{f_ix}, firstLine(f_ix), regFail(f_ix)]= alignAsync(workerTable, params, f_ix, DMD_ix);
        end
    end
    trialTable.motion_correction.registration_failed(DMD_ix,:) = regFail;
    if params.isReVolt && ~isfield(trialTable.motion_correction, 'first_line_original')
        trialTable.motion_correction.first_line_original = trialTable.slap2_info.first_line;
        tOffset = firstLine - trialTable.motion_correction.first_line_original(DMD_ix,:);
        trialTable.slap2_info.first_line = trialTable.slap2_info.first_line + tOffset;
    end
    trialTable.motion_correction.fn_reg_ds(DMD_ix,:) = fnRegDS;
    trialTable.motion_correction.fn_adata(DMD_ix,:) = fnAdata;
end
%during alignment of some data we discard initial frames

params.endTime = char(datetime('now','TimeZone','local','Format','yyyy-MM-dd''T''HH:mm:ss.SSSZZZZZ'));

trialTable.motion_correction.align_params = params;
saveStructToH5(trialTable, [trialtabledr filesep fn]);

disp('done multiRoiRegistration.')
end

function [fnwrite, fnAdata, firstLine, registrationFailed] = alignAsync(trialTable, params, f_ix, DMD_ix)
if params.includeIntegrationROIs
    spTypeFlag = 0; %use all superpixel types
else
    spTypeFlag = 1; %use only raster superpixels
end

mocosavedr = fullfile(trialTable.savedr,'motion_correction');
fn = trialTable.filename{DMD_ix,f_ix};
fnW = ['E' int2str(trialTable.epoch(DMD_ix,f_ix)) 'T' int2str(f_ix) 'DMD' int2str(DMD_ix)];
firstLine = trialTable.slap2_info.first_line(DMD_ix,f_ix);
lastLine = trialTable.slap2_info.last_line(DMD_ix, f_ix);
aData = params;

% Exact MEX/MATLAB dispatch options. The initial-template stage is always
% exhaustive; adaptive search, if explicitly enabled, is restricted to the
% main frame-registration loop.
xcorrOptsExact = makeWeightedXcorrOptions(aData,false);
xcorrOptsMain = makeWeightedXcorrOptions(aData,aData.useAdaptiveWeightedXcorr);
initXcorrStats = newXcorrStats();
mainXcorrStats = newXcorrStats();

disp(['Aligning: ' fnW ' of ' [trialTable.datadr filesep fn]])
fnwrite = [fnW '_REGISTERED_DOWNSAMPLED-' int2str(aData.alignHz) 'Hz.tif'];
fnAdata = [fnW '_ALIGNMENTDATA.h5'];
tifPath = fullfile(mocosavedr,fnwrite);
adataPath = fullfile(mocosavedr,fnAdata);
partialAdataPath = [adataPath '.partial'];
registrationFailed = false;
alignWallStart = tic;

% A trial is complete only when BOTH final outputs exist. If only one side
% exists (for example after an HDF5 failure), discard the stale pair before
% recomputing so a restart can never append to/reuse a partial registered TIFF.
if ~params.overwriteExisting && exist(adataPath,'file') && exist(tifPath,'file')
    disp([fnW ' of ' fn ' is already aligned; skipping' newline 'To force realign, set overwriteExisting=TRUE']);
    return
end
cleanupIncompleteOutputPair(adataPath,tifPath,partialAdataPath);

readerReady = false;
lastReaderError = [];
for retry = 1:10
    try
        fullDatPath = fullfile(trialTable.datadr,fn);
        reuseReader = true;
        if isfield(params,'reuseSlap2Reader')
            reuseReader = logical(params.reuseSlap2Reader);
        end
        tReaderSetup = tic;
        [S2data,meta,readerCacheHit] = getCachedRegistrationReader(fullDatPath,reuseReader);
        readerSetupSeconds = toc(tReaderSetup);
        readerReady = true;
        break
    catch ME
        lastReaderError = ME;
    end
end
if ~readerReady
    fprintf(2,'Failed to initialize Slap2DataReader after 10 attempts for %s\n',fullDatPath);
    rethrow(lastReaderError);
end
linerateHz = 1/meta.linePeriod_s;
dt = linerateHz/aData.alignHz;
numChannels = S2data.numChannels;
minSamps = 15; %minimimum number of samples to include in template

if params.isReVolt
    numChannels = 1;
    redChannel =2;
    if isfield(trialTable, 'motion_correction') && isfield(trialTable.motion_correction, 'first_line_original') && ~isnan(trialTable.motion_correction.first_line_original(DMD_ix,f_ix))
        %first line already adjusted
    else
        % load frames until you see the light turn on on channel 2
        f0 = firstLine+1000; minI = []; maxI = [];
        fEnd = round(0.8*firstLine + 0.2*lastLine);
        nSamps = 15;
        span = fEnd-f0;
        while (fEnd-f0)>(0.4*nSamps*dt)
            frames = round(linspace(f0,fEnd,nSamps));
            [redFrames,~] = getImagesWrapper(S2data,redChannel,frames,ceil(dt),1,spTypeFlag);
            redFrames = normalizeImageStack(redFrames,1,nSamps);
            ii = reshape(mean(redFrames,[1 2 3],'omitnan'),1,[]);
            if isempty(minI)
                minI = min(ii, [], 'omitnan');
                maxI = max(ii, [], 'omitnan');
            end
            ixEnd0 = find(ii>(0.2*minI + 0.8*maxI),1,'first');
            if isempty(ixEnd0)
                warning(['isReVolt flag was set but laser on time could not be detected for file:' fnW ' of ' fn '. skipping...'])
                registrationFailed = true;
                return
            end
            ixEnd = min(numel(ii),ixEnd0+1);
            ix00 = find(ii(1:ixEnd)<(0.65*minI + 0.35*maxI),1,'last');
            if isempty(ix00)
                warning(['isReVolt flag was set but laser transition could not be bracketed for file:' fnW ' of ' fn '. skipping...'])
                registrationFailed = true;
                return
            end
            ix0 = max(1,ix00-2);
            if (frames(ixEnd)-frames(ix0)) >=span
                warning('trouble zooming in on time of light turn on.')
                break %stop zooming in
            else
                span = (frames(ixEnd)-frames(ix0));
            end
            fEnd = frames(ixEnd); f0 = frames(ix0);
        end
        iMid = find(ii>(min(ii)+max(ii))/2,1,'first');
        if isempty(iMid) || iMid < 2
            warning(['isReVolt flag was set but laser midpoint could not be interpolated for file:' fnW ' of ' fn '. skipping...'])
            registrationFailed = true;
            return
        end
        firstLine = round(interp1(ii(iMid + [-1 0]),frames(iMid+[-1 0]),(min(ii)+max(ii))/2));
    end
end

%sanity checks
if isprop(S2data, 'hDataFile')
    assert(length(S2data.hDataFile.fastZs)==1); %single plane acquisitions only
    metaZ = S2data.hDataFile.fastZs;
else
    assert(length(S2data.hMultiDataFiles.fastZs)==1)
    metaZ = S2data.hMultiDataFiles.fastZs;
end

%%%%%Make an initial template
%crosscorrelate each initial frame to each other
disp('generating template')
initFrames =   round(linspace(firstLine, lastLine, 42)); initFrames = initFrames(2:end-1);
nInitFrames = length(initFrames);

if nInitFrames==0
    disp(['File ' fnW ' of ' fn ' was very short! Skipping alignment']);
    registrationFailed = true;
    return
end
% Read all initial-template frames/channels in one Slap2DataReader batch.
% The summed registration image is identical to the previous per-frame,
% per-channel getImage loop, but reader/MEX setup is amortized.
tInitialRead = tic;
[Yinit,freshness] = getImagesWrapper(S2data,1:numChannels,initFrames,ceil(dt),1,spTypeFlag);
initialReadSeconds = toc(tInitialRead);
Yinit = normalizeImageStack(Yinit,numChannels,nInitFrames);
freshness = normalizeFreshnessStack(freshness,nInitFrames);
rawImageSize = [size(Yinit,1),size(Yinit,2)];

Y = reshape(sum(Yinit,3),rawImageSize(1),rawImageSize(2),nInitFrames);
clear Yinit

%make data smaller for alignment
firstTrimRow = find(~all(isnan(Y(:,:,1)),2),1,'first');
lastTrimRow = find(~all(isnan(Y(:,:,1)),2),1,'last');
firstTrimCol = find(~all(isnan(Y(:,:,1)),1),1,'first');
lastTrimCol = find(~all(isnan(Y(:,:,1)),1),1,'last');
if isempty(firstTrimRow) || isempty(lastTrimRow) || isempty(firstTrimCol) || isempty(lastTrimCol)
    error('MultiRoiRegistration:NoRasterPixels', ...
        'Initial SLAP2 template frames contain no finite raster pixels for %s.',fnW);
end
trimRows = firstTrimRow:lastTrimRow;
trimCols = firstTrimCol:lastTrimCol;
Y = Y(trimRows, trimCols,:); freshness = freshness(trimRows, trimCols,:);
sz = size(Y);

R = ones(nInitFrames);
motion = zeros(2,nInitFrames,nInitFrames);
tInitialCorr = tic;
for f1 = 1:nInitFrames
    for f2 = (f1+1):nInitFrames
        if aData.useFastWeightedXcorr
            [motion(:,f1,f2),R(f1,f2),~,xcorrInfo] = xcorr2_nans_weighted_dispatch( ...
                Y(:,:,f2),freshness(:,:,f2),Y(:,:,f1),[0;0],3,xcorrOptsExact);
            initXcorrStats = accumulateXcorrStats(initXcorrStats,xcorrInfo);
        else
            [motion(:,f1,f2),R(f1,f2)] = xcorr2_nans_weighted( ...
                Y(:,:,f2),freshness(:,:,f2),Y(:,:,f1),[0;0],3);
            initXcorrStats.dispatchCalls = initXcorrStats.dispatchCalls + 1;
            initXcorrStats.legacyCalls = initXcorrStats.legacyCalls + 1;
        end
        motion(:,f2,f1) = -motion(:,f1,f2);
        R(f2,f1) = R(f1,f2);
    end
end
initialCorrSeconds = toc(tInitialCorr);
[bestR, maxind] = max(median(R));
frameInds = find(R(:,maxind)>=bestR);

assert(aData.maxshift==round(aData.maxshift), 'params.maxshift must be an integer');
[viewR, viewC] = ndgrid((1:(sz(1)+2*aData.maxshift))-aData.maxshift, (1:(sz(2)+2*aData.maxshift))-aData.maxshift); %view matrices for interpolation
tFrames = nan(2*aData.maxshift+sz(1),2*aData.maxshift+sz(2),numel(frameInds));
baseViewC = viewC(1,:);
baseViewR = viewR(:,1);
for fix = 1:numel(frameInds)
    if aData.useFastInterpolation
        tFrames(:,:,fix) = interpFrameTranslationChannels( ...
            Y(:,:,frameInds(fix)),baseViewC,baseViewR, ...
            -motion(2,frameInds(fix),maxind), ...
            -motion(1,frameInds(fix),maxind),freshness(:,:,frameInds(fix)));
    else
        tFrames(:,:,fix) = interpFrame( ...
            Y(:,:,frameInds(fix)), ...
            baseViewC-motion(2,frameInds(fix),maxind), ...
            baseViewR-motion(1,frameInds(fix),maxind), ...
            freshness(:,:,frameInds(fix)));
    end
end
tSum = sum(tFrames,3, 'omitnan');
tN = sum(~isnan(tFrames),3);
template = sqrt(tSum./tN);
template(tN<minSamps) = nan;

% These arrays are no longer needed once the initial template exists.
clear tFrames R motion frameInds

if params.refStackTemplate
    if params.isReVolt
        error('refStack alignment not implemented for reVolt imaging')
    end
    pathKey = ['Path' int2str(DMD_ix)];
    if isfield(trialTable.slap2_info.ref_stack, pathKey)
        dmdRef = trialTable.slap2_info.ref_stack.(pathKey);
    else
        dmdRef = trialTable.slap2_info.ref_stack.(['DMD' int2str(DMD_ix)]);
    end
    fullTemplate = nan(size(dmdRef.IM,[2 1]));
    fullTemplate((min(trimRows)-aData.maxshift):(max(trimRows)+aData.maxshift),(min(trimCols)-aData.maxshift):(max(trimCols)+aData.maxshift)) = template;

    if numel(dmdRef.channels) == 2
        refStack = (dmdRef.IM(:,:,1:2:end) + dmdRef.IM(:,:,2:2:end))/2;
    else
        refStack = dmdRef.IM;
    end
    templateShifts = xcorr2_nans(fullTemplate,refStack(:,:,floor(end/2)+1)',[0;0],aData.maxshift);
    T0 = imtranslate(permute(refStack,[2 1 3]),[templateShifts(2:-1:1),0]);
    T0 = T0((min(trimRows)-aData.maxshift):(max(trimRows)+aData.maxshift),(min(trimCols)-aData.maxshift):(max(trimCols)+aData.maxshift),:);
    disp('template generated from reference stack')
else
    T0 = template;
end

clear Y Yhp;

initR = 0; initC = 0;
DSframes = ceil(firstLine:dt:lastLine);
nDSframes= length(DSframes); %number of downsampled frames

% For spatial downsampling, used to calculate alignment quality. Keep only
% one ~10 s QC chunk in RAM instead of the entire downsampled movie.
dsTimes = 2;
dsSz = floor(size(template)./(2^dsTimes));
qcOffset = 10;
% Preserve the legacy chunk-edge construction exactly. Note that the
% historical code used N edge points and therefore N-1 actual chunks.
nQCEdgePoints = ceil(nDSframes./(aData.alignHz*10));
qcChunkEdges = round(linspace(1,nDSframes+1,nQCEdgePoints));
nQCChunks = max(0,numel(qcChunkEdges)-1);
recNegErr = nan(1,nDSframes);
qcChunkIx = 1;
if nQCChunks>0
    qcChunkStart = qcChunkEdges(1);
    qcChunkEnd = qcChunkEdges(2)-1;
    A_ds_chunk = nan([dsSz,qcChunkEnd-qcChunkStart+1],'single');
else
    qcChunkStart = [];
    qcChunkEnd = [];
    A_ds_chunk = [];
end

motionDSr = nan(1,nDSframes);
motionDSc = nan(1,nDSframes); %matrices to store the inferred motion
if params.refStackTemplate
    motionDSz = nan(1,nDSframes); %matrices to store the inferred motion
end
aErrorDS = nan(1,nDSframes); %alignment error output by dftregistration

%output TIF
pixelscale = 4e4; %PIXEL SIZE IN DOTS PER CM; 250nm

% Per-channel mean over time (aligned frames), numChannels x H x W.
% H,W match interpFrame output / saved TIFF pages: trimmed crop plus maxshift padding.
szOut = [sz(1) + 2*aData.maxshift, sz(2) + 2*aData.maxshift];
sumMeanIM = zeros(numChannels, szOut(1), szOut(2), 'single');
nMeanIM = zeros(numChannels, szOut(1), szOut(2), 'single');

% MEMORY OPTIMIZATION:
% The previous implementation allocated varFacDS as
% H x W x nDSframes x single in RAM (often several GB PER WORKER).
% Create the final HDF5 dataset up front and stream small frame batches to
% disk. The on-disk schema remains /slap2/varFacDS.
staticSave = struct();
staticSave.numChannels = numChannels;
staticSave.frametime = 1/aData.alignHz;
staticSave.alignHz = aData.alignHz;
staticSave.DSframes = DSframes;
staticSave.slap2 = struct();
staticSave.slap2.Z_depths = metaZ;
staticSave.slap2.cropRow = trimRows(1)-aData.maxshift;
staticSave.slap2.cropCol = trimCols(1)-aData.maxshift;
staticSave.slap2.viewC = viewC;
staticSave.slap2.viewR = viewR;
staticSave.slap2.trimRows = trimRows;
staticSave.slap2.trimCols = trimCols;
saveStructToH5(staticSave, partialAdataPath);

% Use spatially tiled chunks so SILo can read localization tiles directly.
% The old H x W x 1 layout was optimal for frame writes but forced SILo to
% copy/rechunk the entire multi-GB dataset before localization.
if ~isfield(aData,'varFacChunkXY') || isempty(aData.varFacChunkXY)
    aData.varFacChunkXY = 128;
end
varFacBufferFrames = min(8,nDSframes);
varChunkXY = max(16,round(aData.varFacChunkXY));
varFacChunkSize = [min(varChunkXY,szOut(1)), ...
                   min(varChunkXY,szOut(2)), ...
                   varFacBufferFrames];
ensureNumericDataset(partialAdataPath, '/slap2/varFacDS', ...
    [szOut(1), szOut(2), nDSframes], 'single', ...
    'ChunkSize', varFacChunkSize);

% PRE-FLIGHT ALL FINAL HDF5 DATASETS BEFORE EXPENSIVE REGISTRATION.
% The previous streaming implementation created these datasets only after
% ~20+ minutes of work, then immediately called h5write. On network storage
% that exposed an HDF5 metadata visibility failure (e.g. /runtime/wall_s was
% reported missing immediately after h5create). Pre-creating and verifying
% every destination here makes such failures immediate and leaves only
% writes to already-existing datasets during finalization.
precreateFinalAlignmentDatasets(partialAdataPath,numChannels,szOut,nDSframes,params.refStackTemplate);

% Open the TIFF only after the HDF5 schema preflight succeeds, so a schema
% error cannot leave an open/half-created TIFF behind.
fTIF = Fast_BigTiff_Write(tifPath,pixelscale,0);

varFacBuffer = nan(szOut(1), szOut(2), varFacBufferFrames, 'single');
varFacBufferCount = 0;
varFacBufferStart = 1;

% Registration read batching. The adaptive template remains strictly
% sequential; only raw image reconstruction is batched.
requestedRegistrationBlockFrames = max(1,round(aData.registrationBlockFrames));
registrationBlockFrames = requestedRegistrationBlockFrames;
if isfinite(aData.registrationBlockMemoryGB) && aData.registrationBlockMemoryGB>0
    rawPxEstimate = max(1,double(rawImageSize(1))*double(rawImageSize(2)));
    % Conservative double-precision estimate: channels + freshness + 35% overhead.
    bytesPerFrameEstimate = 8*rawPxEstimate*(numChannels+1)*1.35;
    memoryBytes = double(aData.registrationBlockMemoryGB)*(1024^3);
    maxFramesByMemory = max(1,floor(memoryBytes/bytesPerFrameEstimate));
    registrationBlockFrames = min(registrationBlockFrames,maxFramesByMemory);
end

mainReadSeconds = 0;
mainCorrSeconds = 0;
mainInterpSeconds = 0;
templateUpdateSeconds = 0;
tiffWriteSeconds = 0;
h5WriteSeconds = 0;
qcSeconds = 0;
T0Valid = ~isnan(T0);

disp('Registering:');
try
    for blockStart = 1:registrationBlockFrames:nDSframes
        blockIxs = blockStart:min(nDSframes,blockStart+registrationBlockFrames-1);

        tRead = tic;
        [Yblock,freshBlock] = getImagesWrapper( ...
            S2data,1:numChannels,DSframes(blockIxs),ceil(dt),1,spTypeFlag);
        mainReadSeconds = mainReadSeconds + toc(tRead);

        Yblock = normalizeImageStack(Yblock,numChannels,numel(blockIxs));
        freshBlock = normalizeFreshnessStack(freshBlock,numel(blockIxs));
        Yblock = Yblock(trimRows,trimCols,:,:);
        freshBlock = freshBlock(trimRows,trimCols,:);

        for blockLocalIx = 1:numel(blockIxs)
            DSframeIx = blockIxs(blockLocalIx);
            Mchannels = Yblock(:,:,:,blockLocalIx);
            M1 = Mchannels(:,:,1);
            freshness = freshBlock(:,:,blockLocalIx);
            if numChannels==2
                M2 = Mchannels(:,:,2);
                M = sqrt(M1+M2);
            else
                M = sqrt(M1);
            end

        if ~mod(DSframeIx,1000)
            disp([int2str(DSframeIx) ' of ' int2str(nDSframes)]);
        end

        tCorr = tic;
        if params.refStackTemplate
            T = T0(aData.maxshift-initR + (1:sz(1)),aData.maxshift-initC+(1:sz(2)),:);
            [motOutput,corrCoeff] = xcorr2_nans3d(M,T,[0;0],aData.clipShift);
            mainXcorrStats.dispatchCalls = mainXcorrStats.dispatchCalls + 1;
            mainXcorrStats.otherCalls = mainXcorrStats.otherCalls + 1;
            motionDSz(DSframeIx) = motOutput(3);
        else
            % Compose only the currently-used crop of the static and
            % adaptive templates; avoid copying the full padded template.
            tRows = aData.maxshift-initR + (1:sz(1));
            tCols = aData.maxshift-initC + (1:sz(2));
            T0crop = T0(tRows,tCols);
            adaptiveCrop = template(tRows,tCols);
            T = T0crop;
            t0ValidCrop = T0Valid(tRows,tCols);
            adaptiveValid = ~isnan(adaptiveCrop);
            bothValid = t0ValidCrop & adaptiveValid;
            T(bothValid) = (T0crop(bothValid)+adaptiveCrop(bothValid))/2;
            templateOnly = ~t0ValidCrop & adaptiveValid;
            T(templateOnly) = adaptiveCrop(templateOnly);

            if aData.useFastWeightedXcorr
                [motOutput,corrCoeff,~,xcorrInfo] = xcorr2_nans_weighted_dispatch( ...
                    M,freshness,T,[0;0],aData.clipShift,xcorrOptsMain);
                mainXcorrStats = accumulateXcorrStats(mainXcorrStats,xcorrInfo);
            else
                [motOutput,corrCoeff] = xcorr2_nans_weighted( ...
                    M,freshness,T,[0;0],aData.clipShift);
                mainXcorrStats.dispatchCalls = mainXcorrStats.dispatchCalls + 1;
                mainXcorrStats.legacyCalls = mainXcorrStats.legacyCalls + 1;
            end
        end
        mainCorrSeconds = mainCorrSeconds + toc(tCorr);

        motionDSr(DSframeIx) = initR+motOutput(1);
        motionDSc(DSframeIx) = initC+motOutput(2);
        aErrorDS(DSframeIx) = 1-corrCoeff^2;

        % Compute aligned image(s) and variance factor. For two-channel
        % recordings the fast path builds interpolation indices/freshness
        % weights once and applies them to both channels.
        tInterp = tic;
        if aData.useFastInterpolation
            if numChannels==2
                [Aall,Vframe] = interpFrameTranslationChannels( ...
                    Mchannels,baseViewC,baseViewR, ...
                    motionDSc(DSframeIx),motionDSr(DSframeIx),freshness);
                A1 = Aall(:,:,1);
                A2 = Aall(:,:,2);
            else
                [A1,Vframe] = interpFrameTranslationChannels( ...
                    M1,baseViewC,baseViewR, ...
                    motionDSc(DSframeIx),motionDSr(DSframeIx),freshness);
            end
        else
            [A1,Vframe] = interpFrame(M1, ...
                baseViewC+motionDSc(DSframeIx), ...
                baseViewR+motionDSr(DSframeIx),freshness);
            if numChannels==2
                A2 = interpFrame(M2, ...
                    baseViewC+motionDSc(DSframeIx), ...
                    baseViewR+motionDSr(DSframeIx),freshness);
            end
        end
        mainInterpSeconds = mainInterpSeconds + toc(tInterp);

        varFacBufferCount = varFacBufferCount + 1;
        varFacBuffer(:,:,varFacBufferCount) = single(Vframe);
        if varFacBufferCount == varFacBufferFrames || DSframeIx == nDSframes
            tH5 = tic;
            writeHyperslabRobust(partialAdataPath,'/slap2/varFacDS', ...
                varFacBuffer(:,:,1:varFacBufferCount), ...
                [1,1,varFacBufferStart], ...
                [szOut(1),szOut(2),varFacBufferCount]);
            h5WriteSeconds = h5WriteSeconds + toc(tH5);
            varFacBufferStart = DSframeIx + 1;
            varFacBufferCount = 0;
        end
        clear Vframe

        tTif = tic;
        fTIF.WriteIMG(single(A1));
        if numChannels==2
            fTIF.WriteIMG(single(A2));
            A = A1+A2;
        else
            A = A1;
        end
        tiffWriteSeconds = tiffWriteSeconds + toc(tTif);

        a1 = single(A1); m1 = ~isnan(a1); a1(~m1) = 0;
        sumMeanIM(1,:,:) = sumMeanIM(1,:,:) + reshape(a1, [1, szOut(1), szOut(2)]);
        nMeanIM(1,:,:) = nMeanIM(1,:,:) + reshape(single(m1), [1, szOut(1), szOut(2)]);
        if numChannels==2
            a2 = single(A2); m2 = ~isnan(a2); a2(~m2) = 0;
            sumMeanIM(2,:,:) = sumMeanIM(2,:,:) + reshape(a2, [1, szOut(1), szOut(2)]);
            nMeanIM(2,:,:) = nMeanIM(2,:,:) + reshape(single(m2), [1, szOut(1), szOut(2)]);
        end

        % Downsample in space and compute registration QC one ~10 s chunk
        % at a time instead of retaining A_ds for the full acquisition.
        tQc = tic;
        dsTmp = A;
        for dsIx = 1:dsTimes
            dsTmp = dsTmp(1:2:2*floor(end/2),1:2:2*floor(end/2)) + ...
                    dsTmp(1:2:2*floor(end/2),2:2:2*floor(end/2)) + ...
                    dsTmp(2:2:2*floor(end/2),1:2:2*floor(end/2)) + ...
                    dsTmp(2:2:2*floor(end/2),2:2:2*floor(end/2));
        end
        if nQCChunks>0
            A_ds_chunk(:,:,DSframeIx-qcChunkStart+1) = single(dsTmp);
            if DSframeIx == qcChunkEnd
                recNegErr(qcChunkStart:qcChunkEnd) = ...
                    computeRecNegErrChunk(A_ds_chunk,qcOffset);
                qcChunkIx = qcChunkIx+1;
                if qcChunkIx <= nQCChunks
                    qcChunkStart = qcChunkEdges(qcChunkIx);
                    qcChunkEnd = qcChunkEdges(qcChunkIx+1)-1;
                    A_ds_chunk = nan([dsSz,qcChunkEnd-qcChunkStart+1],'single');
                end
            end
        end
        qcSeconds = qcSeconds + toc(tQc);

        tTemplate = tic;
        sel = ~isnan(A);
        tSum(sel) = tSum(sel)+A(sel);
        tN(sel) = tN(sel)+1;
        ready = sel & tN>=minSamps;
        template(ready) = sqrt(tSum(ready)./tN(ready));
        notReady = sel & tN<minSamps;
        template(notReady) = nan;
        templateUpdateSeconds = templateUpdateSeconds + toc(tTemplate);

        initR = max(-aData.maxshift,min(aData.maxshift,round(motionDSr(DSframeIx))));
        initC = max(-aData.maxshift,min(aData.maxshift,round(motionDSc(DSframeIx))));
        end % frames within reader block
    end % reader blocks
catch ME
    disp(ME);
    aData.registrationFailed = true;
    registrationFailed = true;
    try
        fTIF.close;
    catch
    end
    % Partial outputs must never satisfy the "already aligned" check.
    if exist(partialAdataPath, 'file')
        delete(partialAdataPath);
    end
    if exist(tifPath, 'file')
        delete(tifPath);
    end
    disp(['REGISTRATION ERROR OCCURRED FOR FILE: ' fnW ' of ' fn newline 'YOU MAY NEED TO QC THIS FILE!' newline 'CONTINUING...'])
    return
end

fTIF.close;

meanIM = sumMeanIM ./ max(nMeanIM, single(1));
meanIM(nMeanIM < single(1)) = nan;

coreProfiledSeconds = readerSetupSeconds + initialReadSeconds + initialCorrSeconds + ...
    mainReadSeconds + mainCorrSeconds + mainInterpSeconds + ...
    templateUpdateSeconds + tiffWriteSeconds + h5WriteSeconds + qcSeconds;
coreWallSeconds = toc(alignWallStart);

if std(motionDSc)>1.5 || std(motionDSr)>1.5
    aData.registrationFailed = true;
    warning(['Too much motion in file ' fnW ' but will still save file']);
else
    aData.registrationFailed = false;
end

%save alignment metadata
aData.numChannels = numChannels;
aData.frametime = 1/aData.alignHz;
aData.DSframes = DSframes;
aData.motionDSc = motionDSc;
aData.motionDSr = motionDSr;
% varFacDS was streamed directly to /slap2/varFacDS in partialAdataPath.
if params.refStackTemplate
    aData.motionDSz = motionDSz;
end
aData.aError = aErrorDS;
aData.Z_depths = metaZ;
%aData.aRankCorrDS = aRankCorrDS;
aData.recNegErr = recNegErr;
aData.cropRow = trimRows(1)-aData.maxshift; %offset to add to ROIs to index into original recording
aData.cropCol = trimCols(1)-aData.maxshift; %offset to add to ROIs to index into original recording

disp('Getting online motion correction offsets')
tOnlineMotion = tic;
if isprop(S2data, 'hDataFile')
    [aData.onlineXshift, aData.onlineYshift, aData.onlineZshift] = getOnlineMotion(S2data.hDataFile, DSframes);
else
    [aData.onlineXshift, aData.onlineYshift, aData.onlineZshift] = getOnlineMotion(S2data.hMultiDataFiles, DSframes);
end
onlineMotionSeconds = toc(tOnlineMotion);
%CONVERTING DATAFILE IMAGES INTO THE SAVED TIFF IMAGE SPACE:
aData.trimRows = trimRows; %used to remap images from the datafile into the space of the saved tiffs
aData.trimCols = trimCols;%used to remap images from the datafile into the space of the saved tiffs
aData.viewC = viewC;%used to remap images from the datafile into the space of the saved tiffs
aData.viewR = viewR;%used to remap images from the datafile into the space of the saved tiffs
%EXAMPLE CODE
%  Y = S2data.getImage(channel, lineTime, deltaTime, zPos);
%  Ytrimmed = Y(aData.trimRows, aData.trimCols);
%  sz = size(Ytrimmed);
%  Yshifted = interp2(1:sz(2), 1:sz(1), Ytrimmed,aData.viewC+motionC, aData.viewR+motionR, 'linear', nan);

registrationFailed = aData.registrationFailed;

% Complete the partially-written alignment H5. All destination datasets
% were pre-created before registration, so finalization performs writes only.
tFinalH5 = tic;
writeNumericDatasetRobust(partialAdataPath, '/meanIM', meanIM);
writeNumericDatasetRobust(partialAdataPath, '/motionDSc', aData.motionDSc);
writeNumericDatasetRobust(partialAdataPath, '/motionDSr', aData.motionDSr);
writeNumericDatasetRobust(partialAdataPath, '/recNegErr', aData.recNegErr);
writeNumericDatasetRobust(partialAdataPath, '/registrationFailed', aData.registrationFailed);
totalXcorrStats = combineXcorrStats(initXcorrStats,mainXcorrStats);
writeNumericDatasetRobust(partialAdataPath, '/runtime/xcorrInitDispatchCalls', initXcorrStats.dispatchCalls);
writeNumericDatasetRobust(partialAdataPath, '/runtime/xcorrMainDispatchCalls', mainXcorrStats.dispatchCalls);
writeNumericDatasetRobust(partialAdataPath, '/runtime/xcorrMexAttempts', totalXcorrStats.mexAttempts);
writeNumericDatasetRobust(partialAdataPath, '/runtime/xcorrMexSuccesses', totalXcorrStats.mexSuccesses);
writeNumericDatasetRobust(partialAdataPath, '/runtime/xcorrMatlabFastCalls', totalXcorrStats.matlabFastCalls);
writeNumericDatasetRobust(partialAdataPath, '/runtime/xcorrMexFallbacks', totalXcorrStats.mexFallbacks);
writeNumericDatasetRobust(partialAdataPath, '/runtime/xcorrLegacyCalls', totalXcorrStats.legacyCalls);
writeNumericDatasetRobust(partialAdataPath, '/runtime/xcorrOtherCalls', totalXcorrStats.otherCalls);
writeNumericDatasetRobust(partialAdataPath, '/runtime/xcorrMainMexSuccesses', mainXcorrStats.mexSuccesses);
writeNumericDatasetRobust(partialAdataPath, '/runtime/xcorrMainMatlabFastCalls', mainXcorrStats.matlabFastCalls);

writeNumericDatasetRobust(partialAdataPath, '/runtime/readerSetup_s', readerSetupSeconds);
writeNumericDatasetRobust(partialAdataPath, '/runtime/initialRead_s', initialReadSeconds);
writeNumericDatasetRobust(partialAdataPath, '/runtime/initialCorrelation_s', initialCorrSeconds);
writeNumericDatasetRobust(partialAdataPath, '/runtime/getImages_s', mainReadSeconds);
writeNumericDatasetRobust(partialAdataPath, '/runtime/correlation_s', mainCorrSeconds);
writeNumericDatasetRobust(partialAdataPath, '/runtime/interpolation_s', mainInterpSeconds);
writeNumericDatasetRobust(partialAdataPath, '/runtime/templateUpdate_s', templateUpdateSeconds);
writeNumericDatasetRobust(partialAdataPath, '/runtime/tiffWrite_s', tiffWriteSeconds);
writeNumericDatasetRobust(partialAdataPath, '/runtime/h5Write_s', h5WriteSeconds);
writeNumericDatasetRobust(partialAdataPath, '/runtime/qc_s', qcSeconds);
writeNumericDatasetRobust(partialAdataPath, '/runtime/coreWall_s', coreWallSeconds);
writeNumericDatasetRobust(partialAdataPath, '/runtime/onlineMotion_s', onlineMotionSeconds);
writeNumericDatasetRobust(partialAdataPath, '/runtime/registrationBlockFrames', registrationBlockFrames);
writeNumericDatasetRobust(partialAdataPath, '/slap2/onlineMotionXshift', aData.onlineXshift);
writeNumericDatasetRobust(partialAdataPath, '/slap2/onlineMotionYshift', aData.onlineYshift);
writeNumericDatasetRobust(partialAdataPath, '/slap2/onlineMotionZshift', aData.onlineZshift);
if isfield(aData, 'motionDSz')
    writeNumericDatasetRobust(partialAdataPath, '/motionDSz', aData.motionDSz);
end
finalH5Seconds = toc(tFinalH5);
profiledTotalSeconds = coreProfiledSeconds + onlineMotionSeconds + finalH5Seconds;
writeNumericDatasetRobust(partialAdataPath, '/runtime/finalH5_s', finalH5Seconds);
writeNumericDatasetRobust(partialAdataPath, '/runtime/profiledTotal_s', profiledTotalSeconds);
prePublishWallSeconds = toc(alignWallStart);
writeNumericDatasetRobust(partialAdataPath, '/runtime/wall_s', prePublishWallSeconds);

% Publish only a complete metadata file. The registered TIFF was closed
% above; renaming the fully-populated partial H5 is the final commit step.
movefile(partialAdataPath, adataPath, 'f');

fprintf(['Registration timing %s: total %.1f min; core %.1f min; reader %.1f s (%s); ' ...
    'init read %.1f s; init corr %.1f s; getImages %.1f s; xcorr %.1f s; interp %.1f s; ' ...
    'template %.1f s; TIFF %.1f s; streamed H5 %.1f s; QC %.1f s; online motion %.1f s; ' ...
    'final H5 %.1f s; block %d frames; xcorr backend MEX %d/%d success/attempt; ' ...
    'MATLAB-fast %d; MEX fallbacks %d; legacy %d; main MEX %d/%d dispatch\n'], ...
    fnW,prePublishWallSeconds/60,coreWallSeconds/60,readerSetupSeconds, ...
    ternary(readerCacheHit,'cache hit','new reader'),initialReadSeconds,initialCorrSeconds, ...
    mainReadSeconds,mainCorrSeconds,mainInterpSeconds,templateUpdateSeconds, ...
    tiffWriteSeconds,h5WriteSeconds,qcSeconds,onlineMotionSeconds,finalH5Seconds, ...
    registrationBlockFrames,totalXcorrStats.mexSuccesses,totalXcorrStats.mexAttempts, ...
    totalXcorrStats.matlabFastCalls,totalXcorrStats.mexFallbacks,totalXcorrStats.legacyCalls, ...
    mainXcorrStats.mexSuccesses,mainXcorrStats.dispatchCalls);
end

function workerTable = makeWorkerTable(trialTable, params, DMD_ix)
%MAKEWORKERTABLE Construct the minimal trial-table view needed by alignAsync.
% Avoid broadcasting source-extraction state and large reference stacks to
% every process worker when they are not used.

workerTable = struct();
workerTable.savedr = trialTable.savedr;
workerTable.datadr = trialTable.datadr;
workerTable.filename = trialTable.filename;
workerTable.epoch = trialTable.epoch;
workerTable.slap2_info = struct();
workerTable.slap2_info.first_line = trialTable.slap2_info.first_line;
workerTable.slap2_info.last_line = trialTable.slap2_info.last_line;

if params.refStackTemplate
    workerTable.slap2_info.ref_stack = struct();
    pathKey = ['Path' int2str(DMD_ix)];
    dmdKey = ['DMD' int2str(DMD_ix)];
    if isfield(trialTable.slap2_info.ref_stack, pathKey)
        workerTable.slap2_info.ref_stack.(pathKey) = ...
            trialTable.slap2_info.ref_stack.(pathKey);
    elseif isfield(trialTable.slap2_info.ref_stack, dmdKey)
        workerTable.slap2_info.ref_stack.(dmdKey) = ...
            trialTable.slap2_info.ref_stack.(dmdKey);
    else
        error('MultiRoiRegistration:MissingReferenceStack', ...
            'No reference stack found for DMD %d.', DMD_ix);
    end
end

if params.isReVolt && isfield(trialTable, 'motion_correction') && ...
        isfield(trialTable.motion_correction, 'first_line_original')
    workerTable.motion_correction.first_line_original = ...
        trialTable.motion_correction.first_line_original;
end
end


function cleanupIncompleteOutputPair(adataPath,tifPath,partialAdataPath)
%CLEANUPINCOMPLETEOUTPUTPAIR Remove outputs unless the final H5/TIFF pair is complete.
paths = {partialAdataPath,adataPath,tifPath};
for ix = 1:numel(paths)
    if exist(paths{ix},'file')
        delete(paths{ix});
    end
end
end


function precreateFinalAlignmentDatasets(filename,numChannels,szOut,nDSframes,hasMotionZ)
%PREFINALALIGNMENTDATASETS Fail fast if the final H5 schema cannot be created.
ensureNumericDataset(filename,'/meanIM',[numChannels,szOut(1),szOut(2)],'single');
ensureNumericDataset(filename,'/motionDSc',[1,nDSframes],'double');
ensureNumericDataset(filename,'/motionDSr',[1,nDSframes],'double');
ensureNumericDataset(filename,'/recNegErr',[1,nDSframes],'double');
ensureNumericDataset(filename,'/registrationFailed',[1,1],'int8');

runtimeScalars = { ...
    '/runtime/xcorrInitDispatchCalls', ...
    '/runtime/xcorrMainDispatchCalls', ...
    '/runtime/xcorrMexAttempts', ...
    '/runtime/xcorrMexSuccesses', ...
    '/runtime/xcorrMatlabFastCalls', ...
    '/runtime/xcorrMexFallbacks', ...
    '/runtime/xcorrLegacyCalls', ...
    '/runtime/xcorrOtherCalls', ...
    '/runtime/xcorrMainMexSuccesses', ...
    '/runtime/xcorrMainMatlabFastCalls', ...
    '/runtime/readerSetup_s', ...
    '/runtime/initialRead_s', ...
    '/runtime/initialCorrelation_s', ...
    '/runtime/getImages_s', ...
    '/runtime/correlation_s', ...
    '/runtime/interpolation_s', ...
    '/runtime/templateUpdate_s', ...
    '/runtime/tiffWrite_s', ...
    '/runtime/h5Write_s', ...
    '/runtime/qc_s', ...
    '/runtime/coreWall_s', ...
    '/runtime/onlineMotion_s', ...
    '/runtime/finalH5_s', ...
    '/runtime/profiledTotal_s', ...
    '/runtime/wall_s', ...
    '/runtime/registrationBlockFrames'};
for ix = 1:numel(runtimeScalars)
    ensureNumericDataset(filename,runtimeScalars{ix},[1,1],'double');
end

% getOnlineMotion deliberately returns column vectors for compatibility.
ensureNumericDataset(filename,'/slap2/onlineMotionXshift',[nDSframes,1],'double');
ensureNumericDataset(filename,'/slap2/onlineMotionYshift',[nDSframes,1],'double');
ensureNumericDataset(filename,'/slap2/onlineMotionZshift',[nDSframes,1],'double');
if hasMotionZ
    ensureNumericDataset(filename,'/motionDSz',[1,nDSframes],'double');
end
end


function ensureNumericDataset(filename,path,sz,dtype,varargin)
%ENSURENUMERICDATASET Create and verify a numeric H5 dataset with retries.
% High-level HDF5 metadata operations can be briefly inconsistent on SMB
% storage. Verify visibility here, before registration, rather than failing
% after the expensive computation has completed. Additional name/value pairs
% (for example 'ChunkSize') are forwarded to H5CREATE.
maxCreateAttempts = 4;
maxVisibilityChecks = 10;
lastError = [];
lastProbeError = [];

for createAttempt = 1:maxCreateAttempts
    if h5DatasetExistsLocal(filename,path)
        break
    end
    try
        h5create(filename,path,sz,'Datatype',dtype,varargin{:});
    catch ME
        lastError = ME;
    end
    [visible,lastProbeError] = waitForNumericDataset(filename,path,sz,maxVisibilityChecks);
    if visible
        break
    end
    if createAttempt < maxCreateAttempts
        pause(h5BackoffSeconds(createAttempt));
    end
end

[visible,probeError] = waitForNumericDataset(filename,path,sz,maxVisibilityChecks);
if ~visible
    if ~isempty(probeError)
        lastProbeError = probeError;
    end
    if ~isempty(lastError)
        detail = lastError.message;
    elseif ~isempty(lastProbeError)
        % h5create never failed; report WHY the reopen/h5info kept failing
        % instead of the uninformative "never became visible".
        detail = sprintf('h5create succeeded but reopening the dataset failed: %s', ...
            lastProbeError.message);
    else
        detail = 'dataset never became visible after creation';
    end
    error('MultiRoiRegistration:H5PrecreateFailed', ...
        'Failed to create/verify HDF5 dataset %s: %s',path,detail);
end
end


function [tf,lastME] = waitForNumericDataset(filename,path,expectedSize,maxChecks)
tf = false;
lastME = [];
for checkIx = 1:maxChecks
    try
        info = h5info(filename,path);
        actualSize = double(info.Dataspace.Size);
        if ~isequal(actualSize(:)',double(expectedSize(:)'))
            error('MultiRoiRegistration:H5DatasetSizeMismatch', ...
                'Dataset %s has size [%s], expected [%s].',path, ...
                num2str(actualSize),num2str(expectedSize));
        end
        tf = true;
        return
    catch ME
        if strcmp(ME.identifier,'MultiRoiRegistration:H5DatasetSizeMismatch')
            rethrow(ME);
        end
        lastME = ME;
        if checkIx < maxChecks
            pause(h5BackoffSeconds(checkIx));
        end
    end
end
end


function writeHyperslabRobust(filename,path,val,start,count)
%WRITEHYPERSLABROBUST Retry a streamed hyperslab write without recreating data.
% Recreating /slap2/varFacDS after streaming has begun would silently discard
% earlier frames, so a persistent failure is surfaced and the trial is cleaned
% up by the existing registration error handler.
maxAttempts = 8;
lastError = [];
for attempt = 1:maxAttempts
    try
        h5write(filename,path,val,start,count);
        return
    catch ME
        lastError = ME;
        if attempt < maxAttempts
            pause(h5BackoffSeconds(attempt));
        end
    end
end
error('MultiRoiRegistration:H5HyperslabWriteFailed', ...
    'Failed writing %s hyperslab after %d attempts: %s', ...
    path,maxAttempts,lastError.message);
end


function s = h5BackoffSeconds(attempt)
s = min(0.8,0.05*2^(max(0,attempt-1)));
end


function writeNumericDatasetRobust(filename,path,val)
%WRITENUMERICDATASETROBUST Write a pre-created dataset with bounded retries.
if isempty(val)
    return
end
val = gather(val);
if islogical(val)
    val = int8(val);
end

expectedSize = size(val);
maxAttempts = 8;
lastError = [];
for attempt = 1:maxAttempts
    try
        if ~h5DatasetExistsLocal(filename,path)
            % This should have been caught by the preflight, but recover if
            % the backing filesystem lost metadata visibility transiently.
            ensureNumericDataset(filename,path,expectedSize,class(val));
        end
        info = h5info(filename,path);
        actualSize = double(info.Dataspace.Size);
        if ~isequal(actualSize(:)',double(expectedSize(:)'))
            error('MultiRoiRegistration:H5DatasetSizeMismatch', ...
                'Dataset %s has size [%s], value has size [%s].',path, ...
                num2str(actualSize),num2str(expectedSize));
        end
        h5write(filename,path,val);
        return
    catch ME
        lastError = ME;
        if attempt < maxAttempts
            pause(h5BackoffSeconds(attempt));
        end
    end
end
error('MultiRoiRegistration:H5WriteFailed', ...
    'Failed to write HDF5 dataset %s after %d attempts: %s', ...
    path,maxAttempts,lastError.message);
end


function tf = h5DatasetExistsLocal(filename,path)
tf = false;
try
    h5info(filename,path);
    tf = true;
catch
end
end


function meta = loadMetadata(datFilename)
ix = strfind(datFilename, 'DMD'+digitsPattern(1));
metaFilename = [datFilename(1:ix+3) '.meta'];
meta = load(metaFilename, '-mat');
end

function [IM, freshness] = getImageWrapper(S2data, channel, frames, dt, zPlane, spTypeFlag)
if spTypeFlag
    [IM,~,freshness] = S2data.getImage(channel,frames,dt,zPlane,spTypeFlag);
else
    [IM,~,freshness] = S2data.getImage(channel,frames,dt,zPlane);
end
end

function [IM,freshness] = getImagesWrapper(S2data,channels,frames,dt,zPlane,spTypeFlag)
%GETIMAGESWRAPPER Batched equivalent of getImageWrapper.
if spTypeFlag
    [IM,freshness] = S2data.getImages(channels,frames,dt,zPlane,spTypeFlag);
else
    [IM,freshness] = S2data.getImages(channels,frames,dt,zPlane);
end
end

function Y = normalizeImageStack(Y,nChannels,nFrames)
%NORMALIZEIMAGESTACK Ensure getImages output is H x W x C x T.
h = size(Y,1);
w = size(Y,2);
nExpected = h*w*nChannels*nFrames;
if numel(Y) ~= nExpected
    error('MultiRoiRegistration:UnexpectedGetImagesShape', ...
        'getImages returned %d elements; expected H*W*%d*%d.', ...
        numel(Y),nChannels,nFrames);
end
Y = reshape(Y,h,w,nChannels,nFrames);
end

function F = normalizeFreshnessStack(F,nFrames)
%NORMALIZEFRESHNESSSTACK Ensure freshness is H x W x T.
h = size(F,1);
w = size(F,2);
if numel(F) ~= h*w*nFrames
    error('MultiRoiRegistration:UnexpectedFreshnessShape', ...
        'Freshness returned %d elements; expected H*W*%d.',numel(F),nFrames);
end
F = reshape(F,h,w,nFrames);
end

function [S2data,meta,cacheHit] = getCachedRegistrationReader(fullDatPath,reuseReader)
%GETCACHEDREGISTRATIONREADER Worker-local reader/metadata cache.
persistent readerMap metadataMap
if isempty(readerMap)
    readerMap = containers.Map('KeyType','char','ValueType','any');
    metadataMap = containers.Map('KeyType','char','ValueType','any');
end

fullDatPath = char(string(fullDatPath));
if ispc
    key = lower(strrep(fullDatPath,'/','\'));
else
    key = fullDatPath;
end

if ~reuseReader
    S2data = slap2.Slap2DataFile(fullDatPath);
    meta = loadMetadata(fullDatPath);
    cacheHit = false;
    return
end

cacheHit = false;
S2data = [];
if exist('slap2.util.getCachedDataFile','file') == 2
    try
        [S2data,cacheHit] = slap2.util.getCachedDataFile(fullDatPath);
    catch
        S2data = [];
    end
end

if isempty(S2data)
    if isKey(readerMap,key)
        S2data = readerMap(key);
        cacheHit = true;
    else
        S2data = slap2.Slap2DataFile(fullDatPath);
        readerMap(key) = S2data;
    end
end

if isKey(metadataMap,key)
    meta = metadataMap(key);
else
    meta = loadMetadata(fullDatPath);
    metadataMap(key) = meta;
end
end

function rec = computeRecNegErrChunk(chunkData,offset)
%COMPUTERECNEGERRCHUNK Exact per-chunk QC calculation from legacy code.
chunkTemplate = median(chunkData,3,'omitmissing');
nanFrac = mean(isnan(chunkData),3);
chunkTemplate(nanFrac>0.2) = nan;
templateGamma = sqrt(max(0,chunkTemplate)+offset);
rec = sqrt(squeeze(mean(max(0,(chunkTemplate-chunkData)./templateGamma).^2, ...
    [1 2],'omitnan') ./ ...
    mean(max(0,(chunkTemplate./templateGamma).^2 + chunkData*0, ...
    'includenan'),[1 2],'omitnan')));
rec = reshape(rec,1,[]);
end

function opts = makeWeightedXcorrOptions(params,useAdaptive)
%MAKEWEIGHTEDXCORROPTIONS Keep per-frame dispatch options compact/constant.
opts = struct();
opts.useMex = params.useMexWeightedXcorr;
opts.validateMexOnFirstUse = params.validateMexWeightedXcorr;
opts.useAdaptive = logical(useAdaptive);
opts.adaptiveRadius = params.adaptiveXcorrRadius;
opts.adaptiveAuditEvery = params.adaptiveXcorrAuditEvery;
opts.adaptiveMinCorrelation = params.adaptiveXcorrMinCorrelation;
end


function stats = newXcorrStats()
%NEWXCORRSTATS Per-trial backend diagnostics; does not affect registration.
stats = struct( ...
    'dispatchCalls',0, ...
    'mexAttempts',0, ...
    'mexSuccesses',0, ...
    'matlabFastCalls',0, ...
    'mexFallbacks',0, ...
    'legacyCalls',0, ...
    'otherCalls',0);
end


function stats = accumulateXcorrStats(stats,info)
%ACCUMULATEXCORRSTATS Add one dispatcher result, including adaptive subcalls.
stats.dispatchCalls = stats.dispatchCalls + 1;
fields = {'mexAttempts','mexSuccesses','matlabFastCalls','mexFallbacks'};
for k = 1:numel(fields)
    f = fields{k};
    if isfield(info,f) && ~isempty(info.(f))
        stats.(f) = stats.(f) + double(info.(f));
    end
end
end


function out = combineXcorrStats(a,b)
%COMBINEXCORRSTATS Combine initial-template and main-loop diagnostics.
out = newXcorrStats();
fields = fieldnames(out);
for k = 1:numel(fields)
    f = fields{k};
    out.(f) = a.(f) + b.(f);
end
end


function out = ternary(condition,trueValue,falseValue)
if condition
    out = trueValue;
else
    out = falseValue;
end
end
