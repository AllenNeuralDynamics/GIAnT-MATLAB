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

%load the trial Table, which sets correspondences between the two DMDs
trialTable = loadStructFromH5([trialtabledr filesep fn]);

mocosavedr = fullfile(trialTable.savedr,'motion_correction');
if ~exist(mocosavedr,'dir')
    mkdir(mocosavedr)
end

%set up parallelization
p = gcp('nocreate'); % If no pool, do not create new one.
if isempty(p)
    poolsize = 1;
else
    poolsize = p.NumWorkers;
end

core_info = evalc('feature(''numcores'');');
core_info = regexp(core_info,'assigned: \d+ logical cores','match');
core_info = core_info{1};
numLogicalCores = str2num(core_info(11:end-14));

nWorkers = min([params.nWorkers, numel(trialTable.filename), numLogicalCores]);
if poolsize<nWorkers || ~strcmpi(class(p), 'parallel.ProcessPool')
    delete(gcp('nocreate'));
    if params.nWorkers<15
        warning('You are using few parallel workers!');
    end
    parpool('processes',nWorkers); %limit the number of workers to avoid running out of RAM %4-30-24, lowering processes again to prevent another error (18 --> 15)
end
nDMDs = size(trialTable.filename,1);
nTrials = size(trialTable.true_trial_ix,2);
trialTable.motion_correction.registration_failed = false(nDMDs, nTrials);

for DMD_ix = 1:nDMDs
    fnRegDS = cell(1,nTrials);
    fnAdata = cell(1,nTrials);
    firstLine = nan(1,nTrials);
    regFail = false(1, nTrials);
    if nWorkers>1
        parfor f_ix = 1:nTrials
            [fnRegDS{f_ix}, fnAdata{f_ix}, firstLine(f_ix), regFail(f_ix)]= alignAsync(trialTable, params, f_ix, DMD_ix);
        end
    else
        for f_ix = 1:nTrials
            [fnRegDS{f_ix}, fnAdata{f_ix}, firstLine(f_ix), regFail(f_ix)]= alignAsync(trialTable, params, f_ix, DMD_ix);
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

disp(['Aligning: ' fnW ' of ' [trialTable.datadr filesep fn]])
fnwrite = [fnW '_REGISTERED_DOWNSAMPLED-' int2str(aData.alignHz) 'Hz.tif'];
if ~params.saveTif
    fnwrite = strrep(fnwrite, '.tif', '.h5');
end
fnAdata = [fnW '_ALIGNMENTDATA.h5'];
registrationFailed = false;

if ~params.overwriteExisting && exist([mocosavedr filesep fnAdata], 'file') && exist([mocosavedr filesep fnwrite], 'file')
    disp([fnW ' of ' fn ' is already aligned; skipping' newline 'To force realign, pass TRUE as second argument']);
    return
end

retry=0;
while retry<10
    try
        S2data = slap2.Slap2DataFile([trialTable.datadr filesep fn]);
        retry = 10;
    catch
        retry = retry+1;
        continue
    end
end
meta = loadMetadata([trialTable.datadr filesep fn]);
linerateHz = 1/meta.linePeriod_s;
dt = linerateHz/aData.alignHz;
numChannels = S2data.numChannels;
minSamps = 15; %minimimum number of samples to include in template

%undo the online motion correction line by line (raster SLAP2 frames only)
comp = setupCompensation(S2data, params, spTypeFlag, fnW);

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
            ii = nan(1,nSamps);
            for fix = 1:nSamps
                ii(fix) = mean(getImageWrapper(S2data, redChannel, frames(fix), ceil(dt), 1, spTypeFlag), 'all', 'omitnan'); %moving image Ch1
            end
            if isempty(minI)
                minI = min(ii, [], 'omitnan');
                maxI = max(ii, [], 'omitnan');
            end
            ixEnd = min(numel(ii), find(ii>(0.2*minI + 0.8*maxI), 1,'first')+1);
            ix0 = max(1,find(ii(1:ixEnd)<(0.65*minI + 0.35*maxI), 1, 'last')-2);
            if (frames(ixEnd)-frames(ix0)) >=span
                warning('trouble zooming in on time of light turn on.')
                break %stop zooming in
            else
                span = (frames(ixEnd)-frames(ix0));
            end
            if isempty(ix0) || isempty(ixEnd)
                warning(['isReVolt flag was set but laser on time could not be detected for file:' fnW ' of ' fn '. skipping...'])
                registrationFailed = true;
                return
            end
            fEnd = frames(ixEnd); f0 = frames(ix0);
        end
        iMid = find(ii>(min(ii)+max(ii))/2, 1, 'first');
        firstLine = round(interp1(ii(iMid + [-1 0]), frames(iMid+[-1 0]), (min(ii)+max(ii))/2));
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
initRef = frameRefOffsets(comp, initFrames, ceil(dt));
for fix = nInitFrames:-1:1
    for cix = numChannels:-1:1
        [Y(:,:,cix,fix), freshness(:,:,fix)] = getImageWrapper(S2data, cix, initFrames(fix), ceil(dt), 1, spTypeFlag, comp, initRef(fix,:));
    end
end

if params.isReVolt
    Y = squeeze(Y(:,:,1,:));
else
    Y = squeeze(sum(Y,3)); %use both channels;
end

%make data smaller for alignment
if comp.enabled
    %compensated frames can place lines up to a few pixels outside the raster;
    %crop to the raster superpixels so that every trial of a path keeps the
    %same output grid (SILo stacks trials by their top-left corner)
    [trimRows, trimCols] = rasterExtent(S2data);
else
    trimRows = find(~all(isnan(Y(:,:,1)),2), 1, 'first'):find(~all(isnan(Y(:,:,1)),2), 1, 'last');
    trimCols = find(~all(isnan(Y(:,:,1)),1), 1, 'first'):find(~all(isnan(Y(:,:,1)),1), 1, 'last');
end
Y = Y(trimRows, trimCols,:); freshness = freshness(trimRows, trimCols,:);
sz = size(Y);

R = ones(nInitFrames);
motion = zeros(2,nInitFrames,nInitFrames);
for f1 = 1:nInitFrames
    for f2 = (f1+1):nInitFrames
        [motion(:,f1,f2), R(f1,f2)] = xcorr2_nans_weighted(Y(:,:,f2), freshness(:,:,f2), Y(:,:,f1), [0 ; 0], 3); 
        motion(:,f2,f1) = -motion(:,f1,f2);
        R(f2,f1) = R(f1,f2);
    end
end
[bestR, maxind] = max(median(R));
frameInds = find(R(:,maxind)>=bestR);

assert(aData.maxshift==round(aData.maxshift), 'params.maxshift must be an integer');
[viewR, viewC] = ndgrid((1:(sz(1)+2*aData.maxshift))-aData.maxshift, (1:(sz(2)+2*aData.maxshift))-aData.maxshift); %view matrices for interpolation
tFrames = nan(2*aData.maxshift+sz(1), 2*aData.maxshift+sz(2), numel(frameInds));
for fix = 1:numel(frameInds)
    tFrames(:,:,fix) = interpFrame(Y(:,:,frameInds(fix)), viewC(1,:)-motion(2,frameInds(fix), maxind), viewR(:,1)-motion(1,frameInds(fix), maxind), freshness(:,:,frameInds(fix)));
end
tSum = sum(tFrames,3, 'omitnan');
tN = sum(~isnan(tFrames),3);
template = sqrt(tSum./tN);
template(tN<minSamps) = nan;

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
DSref = frameRefOffsets(comp, DSframes, ceil(dt)); %each frame's reference offset when compensating

%for spatial downsampling, used to calculate alignment quality and for quicker visualization
dsTimes = 2;
dsSz = floor(size(template)./(2^dsTimes));
A_ds = nan([dsSz nDSframes],'single');

motionDSr = nan(1,nDSframes);
motionDSc = nan(1,nDSframes); %matrices to store the inferred motion
if params.refStackTemplate
    motionDSz = nan(1,nDSframes); %matrices to store the inferred motion
end
aErrorDS = nan(1,nDSframes); %alignment error output by dftregistration

%output TIF
pixelscale = 4e4; %PIXEL SIZE IN DOTS PER CM; 250nm
if params.saveTif
    fTIF = Fast_BigTiff_Write([mocosavedr filesep fnwrite],pixelscale,0);
end

V1 = nan(size(viewC,1),size(viewC,2),nDSframes,'single'); %variance factor; multiply the image value by this to get variance
%T = T0(aData.maxshift + (1:sz(1)), aData.maxshift+(1:sz(2)));
disp('Registering:');
% Per-channel mean over time (aligned frames), numChannels x H x W.
% H,W match interpFrame output / saved TIFF pages: trimmed crop plus maxshift padding.
szOut = [sz(1) + 2*aData.maxshift, sz(2) + 2*aData.maxshift];
if ~params.saveTif
    h5fn = [mocosavedr filesep fnwrite];
    if exist(h5fn, 'file')
        delete(h5fn);
    end
    h5sz = [szOut nDSframes*numChannels];
    h5create(h5fn, '/data', h5sz, 'Datatype', 'single', 'Deflate', 4, 'ChunkSize', [h5sz(1) h5sz(2) 1]);
end
sumMeanIM = zeros(numChannels, szOut(1), szOut(2), 'single');
nMeanIM = zeros(numChannels, szOut(1), szOut(2), 'single');
try
    for DSframeIx = 1:nDSframes
        [M1, freshness] = getImageWrapper(S2data, 1, DSframes(DSframeIx), ceil(dt), 1, spTypeFlag, comp, DSref(DSframeIx,:)); %moving image Ch1
        M1 = M1(trimRows, trimCols);
        freshness = freshness(trimRows, trimCols); %image freshness, the effective number of samples that contribute to the measurement
        if numChannels==2
            M2 =  getImageWrapper(S2data, 2, DSframes(DSframeIx), ceil(dt), 1, spTypeFlag, comp, DSref(DSframeIx,:)); %moving image Ch2
            M2 = M2(trimRows, trimCols);
            M = sqrt(M1+M2);
        else
            M = sqrt(M1);
        end

        if ~mod(DSframeIx, 1000)
            disp([int2str(DSframeIx) ' of ' int2str(nDSframes)]);
        end

        if params.refStackTemplate
            T = T0(aData.maxshift-initR + (1:sz(1)), aData.maxshift-initC+(1:sz(2)),:);
            [motOutput, corrCoeff] = xcorr2_nans3d(M, T, [0 ; 0], aData.clipShift);
            motionDSz(DSframeIx) = motOutput(3);
        else
            Ttmp = mean(cat(3, T0,template),3, 'omitnan');
            T = Ttmp(aData.maxshift-initR + (1:sz(1)), aData.maxshift-initC+(1:sz(2)));
            %[motOutput, corrCoeff] = xcorr2_nans(M, T, [0 ; 0], aData.clipShift);
            [motOutput, corrCoeff] = xcorr2_nans_weighted(M, freshness, T, [0 ; 0], max(abs([initC initR aData.clipShift])));
        end

        motionDSr(DSframeIx) = initR+motOutput(1);
        motionDSc(DSframeIx) = initC+motOutput(2);
        aErrorDS(DSframeIx) = 1-corrCoeff^2;

        %compute aligned image and variance factor
        [A1,V1(:,:,DSframeIx)] = interpFrame(M1, viewC(1,:)+motionDSc(DSframeIx), viewR(:,1)+motionDSr(DSframeIx), freshness);

        if params.saveTif
            fTIF.WriteIMG(single(A1));
        else
            h5write(h5fn, '/data', single(A1), [1 1 (DSframeIx-1)*numChannels+1], [szOut 1]);
        end
        if numChannels==2
            A2 = interpFrame(M2, viewC(1,:)+motionDSc(DSframeIx), viewR(:,1)+motionDSr(DSframeIx), freshness);%interp2(1:sz(2), 1:sz(1), M2,viewC+motionDSc(DSframeIx), viewR+motionDSr(DSframeIx), 'linear', nan);
            if params.saveTif
                fTIF.WriteIMG(single(A2));
            else
                h5write(h5fn, '/data', single(A2), [1 1 (DSframeIx-1)*numChannels+2], [szOut 1]);
            end
            A =A1+A2;
        else
            A = A1;
        end

        a1 = single(A1); m1 = ~isnan(a1); a1(~m1) = 0;
        sumMeanIM(1,:,:) = sumMeanIM(1,:,:) + reshape(a1, [1, szOut(1), szOut(2)]);
        nMeanIM(1,:,:) = nMeanIM(1,:,:) + reshape(single(m1), [1, szOut(1), szOut(2)]);
        if numChannels==2
            a2 = single(A2); m2 = ~isnan(a2); a2(~m2) = 0;
            sumMeanIM(2,:,:) = sumMeanIM(2,:,:) + reshape(a2, [1, szOut(1), szOut(2)]);
            nMeanIM(2,:,:) = nMeanIM(2,:,:) + reshape(single(m2), [1, szOut(1), szOut(2)]);
        end

        %downsample in space
        dsTmp = A;
        for dsIx = 1:dsTimes
            dsTmp = dsTmp(1:2:2*floor(end/2), 1:2:2*floor(end/2)) + dsTmp(1:2:2*floor(end/2), 2:2:2*floor(end/2)) + dsTmp(2:2:2*floor(end/2), 1:2:2*floor(end/2)) + dsTmp(2:2:2*floor(end/2), 2:2:2*floor(end/2));
        end
        A_ds(:,:,DSframeIx) = dsTmp;

        sel = ~isnan(A);
        tSum(sel) = tSum(sel)+A(sel);
        % tSum = tSum*(1-aData.alpha);
        % tN = tN*(1-aData.alpha);
        tN(sel) = tN(sel)+1;
        template(sel) = sqrt(tSum(sel)./tN(sel));
        template(tN<minSamps) = nan;

        initR = max(-aData.maxshift,min(aData.maxshift, round(motionDSr(DSframeIx))));
        initC = max(-aData.maxshift,min(aData.maxshift,round(motionDSc(DSframeIx))));
    end
catch ME
    disp(ME);
    aData.registrationFailed = true;
    registrationFailed = true;
    try
        fTIF.close;
    catch
    end
    disp(['REGISTRATION ERROR OCCURRED FOR FILE: ' fnW ' of ' fn newline 'YOU MAY NEED TO QC THIS FILE!' newline 'CONTINUING...'])
    return
end

if params.saveTif
    fTIF.close;
end

meanIM = sumMeanIM ./ max(nMeanIM, single(1));
meanIM(nMeanIM < single(1)) = nan;

%Compute alignment error
offset = 10;
recNegErr = nan(1,nDSframes);
nChunks = ceil(nDSframes./(aData.alignHz*10)); %align 10 seconds at a time
chunkEdges = round(linspace(1, nDSframes+1, nChunks+1)); %nChunks chunks of at most 10 s
for chunkIx = 1:length(chunkEdges)-1
    t_ixs = chunkEdges(chunkIx):chunkEdges(chunkIx+1)-1;
    chunkData = A_ds(:,:,t_ixs);
    template = median(chunkData, 3, 'omitmissing');
    nanFrac = mean(isnan(chunkData), 3);
    template(nanFrac>0.2) = nan;
    template_gamma = sqrt(max(0,template)+offset);
    % chunkData*0: broadcast 2D template to 3D without repmat; nan where frame is nan
    recNegErr(1,t_ixs) = sqrt(squeeze(mean(max(0, (template-chunkData)./template_gamma).^2, [1 2], 'omitnan')./mean(max(0, (template./template_gamma).^2 + chunkData * 0, 'includenan'), [1 2], 'omitnan')));
end

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
aData.varFacDS = V1; %variance factor; multiply pixel intensity by this to get a number proportional to the pixel's variance
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
if isprop(S2data, 'hDataFile')
    [aData.onlineXshift, aData.onlineYshift, aData.onlineZshift] = getOnlineMotion(S2data.hDataFile, DSframes);
else
    [aData.onlineXshift, aData.onlineYshift, aData.onlineZshift] = getOnlineMotion(S2data.hMultiDataFiles, DSframes);
end
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

%build the on-disk alignment struct, following the layout documented in README.md
toSave = struct();
toSave.numChannels = aData.numChannels;
toSave.meanIM = meanIM;
toSave.frametime = aData.frametime;
toSave.alignHz = aData.alignHz;
toSave.DSframes = aData.DSframes;
toSave.motionDSc = aData.motionDSc;
toSave.motionDSr = aData.motionDSr;
toSave.recNegErr = aData.recNegErr;
toSave.registrationFailed = aData.registrationFailed;
toSave.slap2 = struct();
toSave.slap2.varFacDS = aData.varFacDS;
toSave.slap2.Z_depths = aData.Z_depths;
toSave.slap2.cropRow = aData.cropRow;
toSave.slap2.cropCol = aData.cropCol;
toSave.slap2.viewC = aData.viewC;
toSave.slap2.viewR = aData.viewR;
toSave.slap2.trimRows = aData.trimRows;
toSave.slap2.trimCols = aData.trimCols;
toSave.slap2.onlineMotionXshift = aData.onlineXshift;
toSave.slap2.onlineMotionYshift = aData.onlineYshift;
toSave.slap2.onlineMotionZshift = aData.onlineZshift;
toSave.slap2.onlineCompensation = compensationRecord(comp, DSref);
if isfield(aData, 'motionDSz')
    toSave.motionDSz = aData.motionDSz;
end
saveStructToH5(toSave, [mocosavedr filesep fnAdata]);
end

function meta = loadMetadata(datFilename)
ix = strfind(datFilename, 'DMD'+digitsPattern(1));
metaFilename = [datFilename(1:ix+3) '.meta'];
meta = load(metaFilename, '-mat');
end

function [IM, freshness] = getImageWrapper(S2data, channel, frames, dt, zPlane, spTypeFlag, comp, refOffset)
%getImage's frame: raster superpixels only (spTypeFlag 1) or all (0). Raster
%frames of SLAP2 recordings come from slap2Frames, which forms the same image
%several times faster than getImage. With comp enabled (see
%setupCompensation), the frame's lines are first moved by their online offsets
%relative to refOffset, the frame's own offset (frameRefOffsets): the frame
%stays in the DMD coordinates of the uncompensated frame, but lines on another
%offset no longer blur it.
if nargin < 7, comp = struct('enabled', false); end
if comp.enabled
    [IM, freshness] = slap2Frames(S2data, channel, frames, dt, ZIdx = zPlane, ...
        CompensateOnline = true, Offsets = comp.offsets, RefOffset = refOffset, Mapping = comp.mapping);
    return
end
if spTypeFlag == 1 && hasLineData(S2data)
    [IM, freshness] = slap2Frames(S2data, channel, frames, dt, ZIdx = zPlane);
    return
end
if spTypeFlag
    [IM,~,freshness] = S2data.getImage(channel, frames, dt, zPlane, spTypeFlag);
else
    [IM,~,freshness] = S2data.getImage(channel, frames, dt, zPlane); %for backward compatibility
end
end

function comp = setupCompensation(S2data, params, spTypeFlag, fnW)
%whether and how to undo the online motion correction (params.compensateOnline):
%only for raster frames of SLAP2 recordings, whose line headers hold the offsets
comp = struct('enabled', false);
if ~isfield(params, 'compensateOnline') || ~params.compensateOnline
    return
end
if spTypeFlag ~= 1 || ~hasLineData(S2data)
    warning('MultiRoiRegistration:noCompensation', ...
        '%s: the online motion correction can only be undone for raster superpixels of SLAP2 recordings; registering without it.', fnW);
    return
end
comp.enabled = true;
comp.offsets = getLineOnlineOffsets(S2data);
comp.refOffset = comp.offsets.medianXY;
comp.mapping = [0 1; 1 0]; %rows follow yOffset_pix, columns xOffset_pix
fprintf('%s: undoing the online motion correction (%d x/y switches)\n', fnW, numel(comp.offsets.xySwitchLines));
end

function [rows, cols] = rasterExtent(S2data)
%the image rows and columns spanned by the raster superpixels of fast-z plane 1
M = S2data.hMultiDataFiles;
NR = double(M.header.dmdPixelsPerRow);
nPx = NR * double(M.header.dmdPixelsPerColumn);
lines = find(M.lineFastZIdxs == 1);
ids = cellfun(@(x) double(x(:)), M.lineSuperPixelIDs(lines), 'UniformOutput', false);
ids = vertcat(ids{:}) + 1;
ids = ids(ids <= nPx);
r = floor((ids - 1) / NR) + 1; c = mod(ids - 1, NR) + 1;
rows = min(r):max(r); cols = min(c):max(c);
end

function tf = hasLineData(S2data)
%whether the reader exposes SLAP2 line data (the real reader always does)
tf = isprop(S2data, 'hMultiDataFiles') && ~isempty(S2data.hMultiDataFiles);
end

function ref = frameRefOffsets(comp, frames, dt)
%each frame's own online offset, nFrames x 2 (empty when not compensating)
ref = zeros(numel(frames), 0);
if comp.enabled
    ref = frameOnlineOffsets(comp.offsets, frames, dt);
end
end

function rec = compensationRecord(comp, DSref)
%the slap2/onlineCompensation group of the alignment file. With compensation,
%each DS frame was formed with its lines moved to the frame's own offset
%frameRefX/Y, and motionDSr/c are relative to that frame's DMD coordinates as
%without compensation. The continuous motion relative to refOffset's frame is
%motionDS + mapping*([frameRefX; frameRefY] - refOffset').
rec = struct('enabled', int8(comp.enabled));
if comp.enabled
    rec.refOffset = comp.refOffset;
    rec.frameRefX = int16(DSref(:,1)');
    rec.frameRefY = int16(DSref(:,2)');
    rec.mapping = comp.mapping;
    rec.yLead = comp.offsets.yLead;
end
end
