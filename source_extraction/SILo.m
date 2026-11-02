function params = SILo(dr_or_pathToTrialTable, paramsIn)
%PARAMETER SETTING
if nargin>1
    if ischar(paramsIn)  % Parse JSON String to Structure
        paramsIn = jsondecode(paramsIn);
    end
    params = setParams('SILo', paramsIn);
else
    params = setParams('SILo');
end
% Fixed values; user/GUI/JSON cannot override these in SILo
params.minBaseline = 0.01;
params.interTrialAlignment = lower(char(params.interTrialAlignment));
if ~any(strcmp(params.interTrialAlignment, {'patchwise', 'affine', 'rigid'}))
    error('SILo:interTrialAlignment', 'params.interTrialAlignment must be ''patchwise'', ''affine'' or ''rigid'', not ''%s''.', params.interTrialAlignment);
end
params.peakMethod = lower(char(params.peakMethod));
if ~any(strcmp(params.peakMethod, {'getactimpeaks', 'wmf'}))
    error('SILo:peakMethod', 'params.peakMethod must be ''getActImPeaks'' or ''wmf'', not ''%s''.', params.peakMethod);
end
if ~nargin
    if ~isInteractiveSession()
        error('SILo:NonInteractive', 'No interactive display: pass the trial_table.h5 path (or its folder) to SILo.');
    end
    [trialTablefn, dr] =  uigetfile('*.h5', 'Select a trial_table file', '*trial_table*.h5' );
else
    if exist(dr_or_pathToTrialTable, 'dir')
        dr = dr_or_pathToTrialTable;
        trialTablefn = 'trial_table.h5';
    else
        [dr trialTablefn ext] = fileparts(dr_or_pathToTrialTable);
        trialTablefn = [trialTablefn ext];
    end
end

params.startTime = char(datetime('now','TimeZone','local','Format','yyyy-MM-dd''T''HH:mm:ss.SSSZZZZZ'));

%confirm that all files exist (also populates source_extraction.fn_raw)
[trialTable, keepTrials] = verifyFiles(trialTablefn, dr);
mocodr = fullfile(trialTable.savedr, 'motion_correction');
nDMDs = size(trialTable.filename,1); % #imaging paths x #trials; non-SLAP2 is treated as 1 path
nTrials = size(trialTable.filename,2);

%parameters that depend only on modality (SLAP2 vs default), hidden from GUI
if ~params.isSLAP2
    params.analyzeHz = nan;
end

disp(['## SUMMARIZING' newline 'Folder:'])
disp(dr)

savedr = fullfile(trialTable.savedr, params.saveSubdir);

if ~exist(savedr, 'dir')
    mkdir(savedr);
end

%call up a GUI for the user to define Soma ROI and regions to exclude
if params.drawUserRois
    %uses existing annotations.h5, else post-acquisition annotations, else draws
    %(drawing only when allowed and someone can respond; otherwise no user ROIs)
    ROIs = loadOrCreateAnnotations(trialTable, keepTrials, mocodr, params.allowAnnotationGUI);
else
    ROIs = [];
end

%PROCESS DATA
for DMDix = nDMDs:-1:1
    firstValidTrial = find(keepTrials(DMDix,:),1,"first");
    if isempty(firstValidTrial)
        warning('SILo:NoValidTrials', ...
            'Skipping imaging path %d: no trials passed file verification.', DMDix);
        exptSummary.Z(DMDix) = nan;
        exptSummary.meanIM{DMDix} = [];
        exptSummary.actIM{DMDix} = [];
        exptSummary.selPix{DMDix} = [];
        exptSummary.sources{DMDix} = struct('R', [], 'C', [], 'V', []);
        exptSummary.userROIs{DMDix} = [];
        exptSummary.peaks{DMDix} = cell(1, nTrials);
        exptSummary.perTrialMeanIMs{DMDix} = [];
        exptSummary.perTrialMeanIMsAligned{DMDix} = [];
        exptSummary.perTrialActIms{DMDix} = [];
        exptSummary.perTrialActIMsAligned{DMDix} = [];
        exptSummary.perTrialAlignmentOffsets{DMDix} = [];
        exptSummary.perTrialAlignmentTforms{DMDix} = [];
        exptSummary.perTrialAlignmentWarp{DMDix} = [];
        continue
    end

    %load some metadata
    fn = trialTable.motion_correction.fn_adata{DMDix,firstValidTrial};
    aData = loadStructFromH5(fullfile(mocodr, fn));
    numChannels = aData.numChannels;
    params.numChannels = numChannels;
    params.alignHz = aData.alignHz;
    if ~params.isSLAP2
        params.analyzeHz = 1/aData.frametime; %analyze conventional recordings at the acquisition framerate
    end
    if isfield(aData, 'slap2') && isfield(aData.slap2, 'Z_depths')
        exptSummary.Z(DMDix) = aData.slap2.Z_depths;
    else
        exptSummary.Z(DMDix) = nan;
        if params.isSLAP2
            warning('Alignment data missing Z_depths, likely out of date!!')
        end
    end
    clear aData

    %set up parallelization
    if params.nWorkers>1
        p = gcp('nocreate');
        if isempty(p)
            poolsize = 0;
        else
            poolsize = p.NumWorkers;
        end
        dd = dir(fullfile(mocodr, trialTable.motion_correction.fn_reg_ds{DMDix, firstValidTrial}));
        try
            fileSize = dd.bytes;
        catch
            error(['Error loading registered tiff:' trialTable.motion_correction.fn_reg_ds{DMDix, firstValidTrial} '\n' 'Are paths in your trial table valid?']);
        end
        if ispc
            userMemInfo = memory;
            memAvailable = userMemInfo.MemAvailableAllArrays;
        elseif ismac
            memAvailable = localMacMemAvailable();
        elseif exist('/proc/meminfo', 'file')
            [~, result] = unix('grep MemAvailable /proc/meminfo | awk ''{print $2}''');
            memKb = str2double(strtrim(result));
            if isnan(memKb)
                warning('SILo:MemProbeFailed', ...
                    'Could not parse MemAvailable from /proc/meminfo; assuming 8 GB.');
                memAvailable = 8 * 1024^3;
            else
                memAvailable = memKb * 1024;  % Convert KB to bytes
            end
        else
            warning('SILo:MemProbeFailed', ...
                'Could not determine available memory; assuming 8 GB.');
            memAvailable = 8 * 1024^3;
        end
        maxWorkers = max(1,min(size(trialTable.filename,2), floor(0.13*memAvailable/fileSize)));
        nWorkers = min(params.nWorkers, maxWorkers);
        
        if poolsize~=nWorkers ||  ~strcmpi(class(p), 'parallel.ProcessPool')
            delete(gcp('nocreate'));
            parpool('processes',nWorkers); %limit the number of workers to avoid running out of RAM
        end
    else
        delete(gcp('nocreate'));
    end

    %Perform Localizations
    disp('Loading data and performing localizations...')
    mIM = cell(1, nTrials); aIM = cell(1,nTrials); alignData = cell(1, nTrials); peaks = cell(1, nTrials); discardFrames = cell(1,nTrials); actAcc = cell(1,nTrials); hIM = cell(1,nTrials);
    fns = trialTable.motion_correction.fn_reg_ds(DMDix, :);
    parfor trialIx = 1:nTrials
        if keepTrials(DMDix,trialIx)
            [~, mIM{trialIx}, aIM{trialIx}, alignData{trialIx}, peaks{trialIx}, discardFrames{trialIx}, actAcc{trialIx}, hIM{trialIx}]= loadAndProcessTrialAsync(mocodr, fns{trialIx}, numChannels, params);
        end
    end
    %Assemble same-sized mean images from different-sized trial means
    szm1 = max(cellfun(@(x)size(x,1),mIM)); szm2 = max(cellfun(@(x)size(x,2), aIM));
    meanIM = nan(szm1,szm2,numChannels, nTrials); activIM = nan(szm1,szm2,1, nTrials);
    halfIM = nan(szm1,szm2,numChannels,2, nTrials); %two independent halves of each trial's mean
    actNumIM = zeros(szm1,szm2,nTrials); actCountIM = zeros(szm1,szm2,nTrials); %activity image sums, 0 where unobserved
    for trialIx = 1:nTrials
        tmp =  mIM{trialIx};
        meanIM(1:size(tmp,1),1:size(tmp,2),:,trialIx) = tmp;
        if ~isempty(hIM{trialIx})
            halfIM(1:size(tmp,1),1:size(tmp,2),:,:,trialIx) = hIM{trialIx};
        end
        tmp =  aIM{trialIx};
        activIM(1:size(tmp,1),1:size(tmp,2),:,trialIx) = tmp;
        if ~isempty(actAcc{trialIx})
            actNumIM(1:size(tmp,1),1:size(tmp,2),trialIx) = actAcc{trialIx}.num;
            actCountIM(1:size(tmp,1),1:size(tmp,2),trialIx) = actAcc{trialIx}.count;
        end
    end

    %Make template
    disp('Making template for aligning across trials...')
    maxshift = 5;
    M = squeeze(sum(meanIM, 3));
    samples = find(keepTrials(DMDix,:)); samples = samples(unique(round(linspace(1,length(samples),20))));
    template = makeConsensusTemplate(M(:,:,samples), maxshift);

    %align all mean images to template: a rigid shift per trial, optionally
    %refined to a patchwise field or an affine map below
    disp('Aligning across trials...')
    corrCoeff = nan(1,nTrials);
    motOutput = nan(2,nTrials);
    Mpad = nan([size(template) size(M,3)]);
    Mpad(maxshift+(1:size(M,1)), maxshift+(1:size(M,2)),:) = M;

    fillval = min(template(:),[], 'omitnan')-1;
    tFFT = fft2(max(template, fillval));
    for trialIx = nTrials:-1:1
        if ~keepTrials(DMDix,trialIx) || all(isnan(activIM(:,:,1,trialIx)), 'all')
            disp(['skipping trial, dmd:' int2str(trialIx) ' ' int2str(DMDix)])
            continue
        end
        disp(['trial: ' int2str(trialIx)])

        output1 = dftregistration_clipped(tFFT, fft2(max(Mpad(:,:,trialIx), fillval)),1,80);
        mot1 = [-output1(3) -output1(4)];
        [motOutput(:,trialIx), corrCoeff(trialIx)] = xcorr2_nans(Mpad(:,:,trialIx), template, round(mot1'), maxshift);
    end
    %transforms from the common frame [row; col; 1] to each trial's image [row; col]
    alignedTrials = find(~isnan(motOutput(1,:)));
    tforms = nan(2,3,nTrials);
    for trialIx = alignedTrials
        tforms(:,:,trialIx) = [1 0 motOutput(1,trialIx); 0 1 motOutput(2,trialIx)];
    end

    %identify outliers in alignment quality to determine valid trials
    corrThresh = min(0.90, median(corrCoeff, 'omitnan')-2*std(corrCoeff, 'omitmissing'));
    actValidPix = squeeze(mean(~isnan(activIM(:,:,1,:)), [1 2]));
    validTrials= find(corrCoeff(:)>corrThresh & actValidPix(:)>mean(actValidPix)/2);

    %trial mean images differ by small deformations that track motion (up to
    %~1 px at the edges of the field, and local bends of a few tenths of a
    %pixel), which a single shift per trial leaves in place; refine the valid
    %trials' alignment groupwise. warp holds every trial's map from the common
    %frame to its image (see trialWarpCoords); trials that are aligned but not
    %valid keep their rigid shift
    warp = struct('affine', tforms, 'nodes', zeros(0,2), 'sigma', nan, 'disp', zeros(2,0,nTrials), 'patch', []);
    if ~strcmpi(params.interTrialAlignment, 'rigid') && numel(validTrials) > 1
        tformsValid = nan(size(tforms));
        tformsValid(:,:,validTrials) = tforms(:,:,validTrials);
        if strcmpi(params.interTrialAlignment, 'patchwise')
            disp('Refining across-trial alignment to patchwise displacement fields...')
            Mhalves = permute(reshape(sum(halfIM, 3), szm1, szm2, 2, nTrials), [1 2 4 3]); %rows x cols x trials x 2, channels summed as in M
            [pw, pwInfo] = alignTrialsPatchwise(M, tformsValid, Mhalves, params.alignPatchPx);
            clear Mhalves
            warp.nodes = pw.nodes; warp.sigma = pw.sigma; warp.patch = pw.patch;
            warp.disp = nan(2, size(pw.nodes,1), nTrials);
            for trialIx = alignedTrials
                warp.disp(:,:,trialIx) = tforms(:,1:2,trialIx)*pw.nodes' + tforms(:,3,trialIx) - pw.nodes';
            end
            warp.disp(:,:,validTrials) = pw.disp(:,:,validTrials);
            warp.affine(:,:,validTrials) = pw.affine(:,:,validTrials);
            disp(['  ' int2str(size(pw.nodes,1)) ' patches of ' int2str(min(pwInfo.nPx)) '-' int2str(max(pwInfo.nPx)) ...
                ' px; median reliability weight ' num2str(median(pwInfo.weight), 2)]);
            disp(['  RMS change from rigid alignment per valid trial (px): ' mat2str(round(pwInfo.rmsDisp(validTrials), 2))]);
        else
            disp('Refining across-trial alignment to affine maps...')
            [tformsValid, affInfo] = alignTrialsAffine(M, tformsValid);
            warp.affine(:,:,validTrials) = tformsValid(:,:,validTrials);
            disp(['  RMS change from rigid alignment per valid trial (px): ' mat2str(round(affInfo.rmsDisp(validTrials), 2))]);
        end
    end
    clear halfIM hIM
    %rigid-equivalent offsets (the displacement at the image centre), used for per-frame motion records
    x0 = [(szm1+1)/2; (szm2+1)/2];
    for trialIx = alignedTrials
        [rC, cC] = trialWarpCoords(warp, trialIx, x0(1), x0(2));
        motOutput(:,trialIx) = [rC; cC] - x0;
    end

    %warp per-trial images into the common frame
    meanAligned = nan(size(meanIM));
    actAligned = nan(size(meanIM,1), size(meanIM,2),1,nTrials);
    for trialIx = alignedTrials
        meanAligned(:,:,:,trialIx) = warpTrialImage(meanIM(:,:,:,trialIx), warp, trialIx, 'cubic');
        actAligned(:,:,1,trialIx) = warpTrialImage(activIM(:,:,1,trialIx), warp, trialIx, 'linear');
    end
    exptSummary.meanIM{DMDix} = mean(meanAligned(:,:,:,validTrials),4, 'omitnan');

    %select sources on the aligned activity image, accumulated over trials in
    %the common frame: sums of local-maximum scores over sums of observations,
    %so pixels observed in only some frames or trials are weighted accordingly
    actNum = zeros(szm1, szm2); actCount = zeros(szm1, szm2); actReg = 0; nScored = 0;
    for trialIx = validTrials(:)'
        num = warpTrialImage(actNumIM(:,:,trialIx), warp, trialIx, 'linear');
        cnt = warpTrialImage(actCountIM(:,:,trialIx), warp, trialIx, 'linear');
        num(isnan(num)) = 0; cnt(isnan(cnt)) = 0;
        actNum = actNum + num;
        actCount = actCount + cnt;
        actReg = actReg + actAcc{trialIx}.reg;
        nScored = nScored + actAcc{trialIx}.nFrames;
    end
    actIM = actNum ./ (actReg + actCount);
    actIM = actIM ./ 10^(floor(log10(max(actIM(:))))-1);
    actIM(actCount < (1-params.nanThresh)*nScored) = nan; %observed in too few of the scored frames
    medIM = nanmedfilt2(actIM, 5.*[1 1]);
    actIM = actIM-medIM; %subtract a local baseline
    exptSummary.actIM{DMDix} = actIM;
    sz = size(actIM);

    %Mask out somata from activity image
    somaMask = false(size(actIM));
    if ~isempty(ROIs)
        for rix = 1:numel(ROIs(DMDix).roiData)
            if contains(upper(ROIs(DMDix).roiData{rix}.Label), 'SOMA')
                tmp = ROIs(DMDix).roiData{rix}.mask;
                somaMask(1:size(tmp,1), 1:size(tmp,2)) = somaMask(1:size(tmp,1), 1:size(tmp,2)) | tmp;
            end
        end
    end

    if strcmp(params.peakMethod, 'wmf')
        %seeds from the trials' aligned activity images rather than the accumulated act_im
        [seedRC, seedInfo] = findSourcesWMF(reshape(actAligned(:,:,1,validTrials), sz(1), sz(2), []), params);
        inSoma = somaMask(sub2ind(sz, round(seedRC(:,1)), round(seedRC(:,2))));
        thetaf = [seedInfo.value(~inSoma), seedRC(~inSoma,:)];
    else
        thetaf = getActImPeaks(actIM,params.peakth,somaMask,params.minPeakDistance);
    end

    sources = struct('R', [], 'C', [], 'V', []);
    totalPix = sum(~isnan(actIM(:)) & ~somaMask(:));
    if totalPix == 0 || isempty(thetaf)
        k = 0;
    else
        sources.R = round(thetaf(:,2));
        sources.C = round(thetaf(:,3));
        sources.V = thetaf(:,1);
        
        % prune source seeds that are the same
        [~, uniqueIx] = unique([sources.R(:) sources.C(:)], 'rows', 'stable');
        sources.R = sources.R(uniqueIx);
        sources.C = sources.C(uniqueIx);
        sources.V = sources.V(uniqueIx);
        k = length(sources.R);
    end

    %select regions near synapses, aligned across movies
    selPix = false([sz(1:2) k]);
    params.selRadius = ceil(2*params.dXY);
    for sourceIx = k:-1:1
        rr = round(sources.R(sourceIx));
        cc = round(sources.C(sourceIx));
        selPix(rr,cc,sourceIx) = true;
        selPix(:,:,sourceIx) = imdilate(selPix(:,:,sourceIx), strel('disk',params.selRadius));
    end
    pxAlwaysValid = mean(isnan(meanAligned(:,:,1,validTrials)),4)<params.nanThresh;
    selPix = selPix & repmat(pxAlwaysValid, 1, 1, k); %exclude poorly measured pixels

    %prune any sources that got clipped by pixel selection
    centerValid = pxAlwaysValid(sub2ind(size(pxAlwaysValid), sources.R, sources.C));
    keepSources = squeeze(sum(selPix, [1 2]))>5 & centerValid(:);
    if k > 0
        sources.R = sources.R(keepSources);
        sources.C = sources.C(keepSources);
        sources.V = sources.V(keepSources);
    end
    selPix = selPix(:,:,keepSources);
    disp(['Number of sources: ' int2str(sum(keepSources))]);
        
    %for each file, load high res data and refine
    params.tau_full=params.tau_s*params.analyzeHz;
    params = setParamsExtractTrial(params);
    
    if isempty(ROIs) || isempty(ROIs(DMDix))
        roiData =[];
    else
        roiData = ROIs(DMDix).roiData;
    end

    if any(keepSources)
        fns = trialTable.source_extraction.fn_raw(DMDix,:);
        if params.isSLAP2
            if isfield(trialTable, 'datadr') && ~isempty(trialTable.datadr)
                trialDataDr = trialTable.datadr;
            else
                trialDataDr = dr;
            end
            fls = trialTable.slap2_info.first_line(DMDix,:);
            els = trialTable.slap2_info.last_line(DMDix,:);
            E = processAllTrials_Async(trialDataDr, fns, fls, els, selPix, sources, discardFrames, alignData, meanAligned, warp, roiData, validTrials, params);
        else % non-SLAP2 (registered movies live under motion_correction)
            fls = cell(1,numel(fns));
            els = cell(1,numel(fns));
            E = processAllTrials_Async(mocodr, fns, fls, els, selPix, sources, discardFrames, alignData, meanAligned, warp, roiData, validTrials, params);
        end
        exptSummary.E(:,DMDix) = E;
    end
    exptSummary.selPix{DMDix} = any(selPix,3);
    exptSummary.sources{DMDix} = sources;
    exptSummary.aData(:,DMDix) = alignData;
    exptSummary.userROIs{DMDix} = roiData;
    exptSummary.peaks{DMDix}= peaks;
    exptSummary.perTrialMeanIMs{DMDix} = meanIM;
    exptSummary.perTrialMeanIMsAligned{DMDix} = meanAligned;
    exptSummary.perTrialActIms{DMDix} = actIM;
    exptSummary.perTrialActIMsAligned{DMDix} = actAligned;
    exptSummary.perTrialAlignmentOffsets{DMDix} = motOutput;
    exptSummary.perTrialAlignmentTforms{DMDix} = warp.affine;
    exptSummary.perTrialAlignmentWarp{DMDix} = warp;

    clear meanAligned meanIM actAligned actNumIM actCountIM E
end

% Shut down the parallel pool explicitly so thread-pool arrays are
% fully materialised into regular MATLAB memory before the HDF5 saves.
delete(gcp('nocreate'));
pause(5);

params.endTime = char(datetime('now','TimeZone','local','Format','yyyy-MM-dd''T''HH:mm:ss.SSSZZZZZ'));

trialTable.source_extraction.analysis_params = params;

exptSummary.params = params;
exptSummary.trialTable = trialTable;
exptSummary.dr = dr;

%save (with retry so a transient HDF5 error does not discard all results)
trySave(@() saveStructToH5(trialTable, [dr filesep trialTablefn]),          'trial_table');
trySave(@() saveExperimentSummaryH5(fullfile(savedr, 'experiment_summary.h5'), exptSummary, trialTable), 'experiment_summary');
trySave(@() savePerTrialSummaryH5(  fullfile(savedr, 'per_trial_summary.h5'),  exptSummary, trialTable), 'per_trial_summary');

%verify the saved file is readable and non-empty
expSumFn = fullfile(savedr, 'experiment_summary.h5');
d = dir(expSumFn);
if isempty(d) || d.bytes == 0
    fprintf(2, 'WARNING: experiment_summary.h5 is missing or empty after save.\n');
else
    fprintf('experiment_summary.h5 written (%d MB)\n', round(d.bytes/1e6));
    for DMDix = 1:numel(exptSummary.sources)
        if isempty(exptSummary.sources{DMDix}) || ~isfield(exptSummary.sources{DMDix}, 'R') ...
                || isempty(exptSummary.sources{DMDix}.R)
            continue
        end
        datasetPath = sprintf('/Path%d/sources/temporal/dF_denoised', DMDix);
        try
            info = h5info(expSumFn, datasetPath);
            fprintf('  %s shape: %s\n', datasetPath, mat2str(info.Dataspace.Size));
        catch
            fprintf('  %s not found.\n', datasetPath);
        end
    end
end


disp('Done SILo')
end

function memAvailable = localMacMemAvailable()
%LOCALMACMEMAVAILABLE Estimate free+inactive memory on macOS via vm_stat.
memAvailable = NaN;
[status, vm] = unix('vm_stat');
if status == 0
    pageSize = 4096;
    tok = regexp(vm, 'page size of (\d+)', 'tokens', 'once');
    if ~isempty(tok)
        pageSize = str2double(tok{1});
    end
    freeTok = regexp(vm, 'Pages free:\s+(\d+)', 'tokens', 'once');
    inactiveTok = regexp(vm, 'Pages inactive:\s+(\d+)', 'tokens', 'once');
    speculativeTok = regexp(vm, 'Pages speculative:\s+(\d+)', 'tokens', 'once');
    if ~isempty(freeTok)
        nPages = str2double(freeTok{1});
        if ~isempty(inactiveTok)
            nPages = nPages + str2double(inactiveTok{1});
        end
        if ~isempty(speculativeTok)
            nPages = nPages + str2double(speculativeTok{1});
        end
        memAvailable = nPages * pageSize;
    end
end
if isnan(memAvailable) || memAvailable <= 0
    [status, result] = unix('sysctl -n hw.memsize');
    memAvailable = str2double(strtrim(result));
    if status ~= 0 || isnan(memAvailable)
        warning('SILo:MemProbeFailed', ...
            'Could not determine available memory on macOS; assuming 8 GB.');
        memAvailable = 8 * 1024^3;
    end
end
end

function trySave(saveFcn, label, maxAttempts)
%TRYSAVE  Call saveFcn up to maxAttempts times, pausing between failures.
% Note: retries only help with transient I/O errors. Failures caused by
% stale thread-pool state in the current process will not recover across
% retries.
if nargin < 3 || isempty(maxAttempts)
    maxAttempts = 3;
end
for attempt = 1:maxAttempts
    try
        saveFcn();
        return;
    catch ME
        fprintf(2, 'WARNING: Attempt %d/%d to save "%s" failed: %s\n', ...
            attempt, maxAttempts, label, ME.message);
        if attempt < maxAttempts
            fprintf(2, 'Waiting 30s before retry...\n');
            pause(30);
        else
            rethrow(ME);
        end
    end
end
end
