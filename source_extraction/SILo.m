function params = SILo(dr_or_pathToTrialTable, paramsIn)
% Selectable source-detection controls are handled locally so existing
% setParams('SILo', ...) behavior remains backward compatible.
requestedSourceDetectionMethod = 'silo';
requestedMaxSynapseDensity = [];
if nargin > 1
    tmpParamsIn = paramsIn;
    if ischar(tmpParamsIn) || (isstring(tmpParamsIn) && isscalar(tmpParamsIn))
        tmpParamsIn = jsondecode(char(tmpParamsIn));
    end
    if isstruct(tmpParamsIn)
        if isfield(tmpParamsIn, 'sourceDetectionMethod')
            requestedSourceDetectionMethod = tmpParamsIn.sourceDetectionMethod;
            tmpParamsIn = rmfield(tmpParamsIn, 'sourceDetectionMethod');
        end
        if isfield(tmpParamsIn, 'maxSynapseDensity')
            requestedMaxSynapseDensity = tmpParamsIn.maxSynapseDensity;
            tmpParamsIn = rmfield(tmpParamsIn, 'maxSynapseDensity');
        end
        paramsIn = tmpParamsIn;
    end
end
%PARAMETER SETTING
if nargin>1
    if ischar(paramsIn)  % Parse JSON String to Structure
        paramsIn = jsondecode(paramsIn);
    end
    params = setParams('SILo', paramsIn);
else
    params = setParams('SILo');
end
params.sourceDetectionMethod = lower(strrep(char(requestedSourceDetectionMethod), '-', '_'));
if any(strcmp(params.sourceDetectionMethod, {'summarizeloco','loco'}))
    params.sourceDetectionMethod = 'summarize_loco';
end
if ~any(strcmp(params.sourceDetectionMethod, {'silo','summarize_loco'}))
    error('SILo:UnknownSourceDetectionMethod', ...
        'sourceDetectionMethod must be "silo" or "summarize_loco".');
end
if ~isempty(requestedMaxSynapseDensity)
    params.maxSynapseDensity = requestedMaxSynapseDensity;
end
% Keep stable copies because setParamsExtractTrial may normalize/drop
% SILo-level parameters that are not used directly by extractTrial.
resolvedSourceDetectionMethod = params.sourceDetectionMethod;
resolvedMaxSynapseDensity = [];
if isfield(params, 'maxSynapseDensity')
    resolvedMaxSynapseDensity = params.maxSynapseDensity;
end
resolvedSavePerTrialSummary = params.savePerTrialSummary;
% Fixed values; user/GUI/JSON cannot override these in SILo
params.minBaseline = 0.01;
runTimer = tic;

% RAM-optimized source-localization settings. These are hidden performance
% controls rather than scientific parameters; callers may override them in
% paramsIn without changing source-localization semantics.
if ~isfield(params, 'localizationTileSize') || isempty(params.localizationTileSize)
    params.localizationTileSize = 96;
end
if ~isfield(params, 'localizationTempDir') || isempty(params.localizationTempDir)
    params.localizationTempDir = tempdir;
end
if ~nargin
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

savedr = fullfile(trialTable.savedr, 'source_extraction');

if ~exist(savedr, 'dir')
    mkdir(savedr);
end

%call up a GUI for the user to define Soma ROI and regions to exclude
if params.drawUserRois
    annotationsDr = fullfile(trialTable.savedr, 'annotations');
    if ~exist(annotationsDr, 'dir')
        mkdir(annotationsDr);
    end
    fnAnnH5 = fullfile(annotationsDr, 'annotations.h5');
    if exist(fnAnnH5, 'file')
        ROIs = loadAnnotationsH5(fnAnnH5);
    else
        drawnDMDs = false(1, nDMDs);
        for DMDix = 1:nDMDs
            firstValidTrial = find(keepTrials(DMDix,:),1,"first");
            if isempty(firstValidTrial)
                warning('SILo:NoValidTrials', ...
                    'Skipping ROI drawing for imaging path %d: no trials passed file verification.', DMDix);
                ROIs(DMDix).dr = mocodr;
                ROIs(DMDix).fn = '';
                ROIs(DMDix).roiData = [];
                continue
            end
            %load image data (prefer meanIM from alignment H5)
            [~, fn, ext] = fileparts(trialTable.motion_correction.fn_reg_ds{DMDix,firstValidTrial});
            adataFn = trialTable.motion_correction.fn_adata{DMDix, firstValidTrial};
            % Do not load /slap2/varFacDS merely to draw annotations.
            % Alignment files can contain a multi-GB variance-factor movie.
            aDataAnnot = loadAlignmentDataLite(fullfile(mocodr, adataFn));
            if isfield(aDataAnnot, 'meanIM') && ~isempty(aDataAnnot.meanIM)
                IM = squeeze(mean(aDataAnnot.meanIM, 1, 'omitnan')); % numChannels x H x W -> H x W
            else
                [IMtif, ~] = ScanImageTiffWrapper(fullfile(mocodr, [fn ext]));
                IM = squeeze(mean(IMtif, [3 4], 'omitnan'));
            end
            hROIs(DMDix) = drawROIs(sqrt(max(0,IM)), savedr, fn);
            ROIs(DMDix).dr = mocodr;
            ROIs(DMDix).fn = fn;
            drawnDMDs(DMDix) = true;
        end
        for DMDix = find(drawnDMDs)
            waitfor(hROIs(DMDix).hF);
            ROIs(DMDix).roiData = hROIs(DMDix).roiData;
        end
        saveAnnotationsH5(fnAnnH5, ROIs); clear hROIs;
    end
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
        exptSummary.rawActIM{DMDix} = [];
        exptSummary.sourceDetectionInfo{DMDix} = struct();
        exptSummary.selPix{DMDix} = [];
        exptSummary.sources{DMDix} = struct('R', [], 'C', [], 'V', []);
        exptSummary.userROIs{DMDix} = [];
        exptSummary.peaks{DMDix} = cell(1, nTrials);
        exptSummary.perTrialMeanIMs{DMDix} = [];
        exptSummary.perTrialMeanIMsAligned{DMDix} = [];
        exptSummary.perTrialActIms{DMDix} = [];
        exptSummary.perTrialActIMsAligned{DMDix} = [];
        exptSummary.perTrialAlignmentOffsets{DMDix} = [];
        continue
    end

    %load some metadata
    fn = trialTable.motion_correction.fn_adata{DMDix,firstValidTrial};
    % Lightweight metadata read: skip /slap2/varFacDS.
    aData = loadAlignmentDataLite(fullfile(mocodr, fn));
    numChannels = aData.numChannels;
    params.numChannels = numChannels;
    if params.activityChannel < 1 || params.activityChannel > numChannels || params.activityChannel ~= round(params.activityChannel)
        error('SILo:InvalidActivityChannel', ...
            'activityChannel=%g is invalid for imaging path %d with %d channel(s).', ...
            params.activityChannel,DMDix,numChannels);
    end
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
    % The RAM-optimized localizer never loads /slap2/varFacDS as a full
    % H x W x T array. The dominant per-worker allocation is therefore the
    % registered TIFF plus bounded tile workspaces. Estimate a conservative
    % worker count from TIFF size and available RAM, and also respect the
    % Processes profile NumWorkers ceiling.
    nWorkers = 1;
    if params.nWorkers>1
        p = gcp('nocreate');
        if isempty(p)
            poolsize = 0;
        else
            poolsize = p.NumWorkers;
        end

        dd = dir(fullfile(mocodr, ...
            trialTable.motion_correction.fn_reg_ds{DMDix, firstValidTrial}));
        if isempty(dd)
            error(['Error loading registered tiff:' ...
                trialTable.motion_correction.fn_reg_ds{DMDix, firstValidTrial} ...
                newline 'Are paths in your trial table valid?']);
        end
        fileSize = double(dd.bytes);

        if ispc
            userMemInfo = memory;
            memAvailable = double(userMemInfo.MemAvailableAllArrays);
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
                memAvailable = memKb * 1024;
            end
        else
            warning('SILo:MemProbeFailed', ...
                'Could not determine available memory; assuming 8 GB.');
            memAvailable = 8 * 1024^3;
        end

        % During TIFF loading, a deinterleaved activity-channel copy can
        % briefly coexist with the full registered movie. Tile-localization
        % then adds bounded work arrays. Reserve ~2.5x the TIFF size plus
        % 1 GiB per process and use no more than 75% of currently available
        % memory for localization workers.
        estimatedWorkerBytes = 2.5*fileSize + 1*1024^3;
        memoryBudget = 0.75*memAvailable;
        maxWorkersByMem = max(1, floor(memoryBudget/estimatedWorkerBytes));

        try
            procCluster = parcluster('Processes');
            maxWorkersByProfile = procCluster.NumWorkers;
        catch
            maxWorkersByProfile = inf;
        end

        nWorkers = min([params.nWorkers, ...
            sum(keepTrials(DMDix,:)), ...
            maxWorkersByMem, ...
            maxWorkersByProfile]);
        nWorkers = max(1, floor(nWorkers));

        fprintf(['Localization workers: %d requested, %d selected ' ...
            '(RAM cap=%d, profile cap=%g)\n'], ...
            params.nWorkers, nWorkers, maxWorkersByMem, maxWorkersByProfile);

        if nWorkers>1
            if poolsize~=nWorkers || ~isa(p, 'parallel.ProcessPool')
                delete(gcp('nocreate'));
                parpool('processes',nWorkers);
            else
                % Reused pool may predate the current code; see
                % refreshPoolWorkerCode.
                refreshPoolWorkerCode(p);
            end
        else
            delete(gcp('nocreate'));
        end
    else
        delete(gcp('nocreate'));
    end

    %Perform Localizations
    disp('Loading data and performing localizations...')
    localizationTimer = tic;
    mIM = cell(1, nTrials); aIM = cell(1,nTrials); alignData = cell(1, nTrials); peaks = cell(1, nTrials); discardFrames = cell(1,nTrials);
    fns = trialTable.motion_correction.fn_reg_ds(DMDix, :);
    if nWorkers>1
        parfor trialIx = 1:nTrials
            if keepTrials(DMDix,trialIx)
                [mIM{trialIx}, aIM{trialIx}, alignData{trialIx}, peaks{trialIx}, discardFrames{trialIx}] = ...
                    loadAndProcessTrialAsync(mocodr, fns{trialIx}, numChannels, params);
            end
        end
    else
        % Avoid MATLAB silently auto-starting a default parallel pool when
        % the memory model has selected a single worker.
        for trialIx = 1:nTrials
            if keepTrials(DMDix,trialIx)
                [mIM{trialIx}, aIM{trialIx}, alignData{trialIx}, peaks{trialIx}, discardFrames{trialIx}] = ...
                    loadAndProcessTrialAsync(mocodr, fns{trialIx}, numChannels, params);
            end
        end
    end
    fprintf('Path %d trial localization: %.1f min\n',DMDix,toc(localizationTimer)/60);

    %Assemble same-sized mean images from different-sized trial means
    szm1 = max(cellfun(@(x)size(x,1),mIM)); szm2 = max(cellfun(@(x)size(x,2), aIM));
    % Visualization/localization images do not require double precision.
    meanIM = nan(szm1,szm2,numChannels,nTrials,'single');
    activIM = nan(szm1,szm2,1,nTrials,'single');
    for trialIx = 1:nTrials
        tmp =  mIM{trialIx};
        meanIM(1:size(tmp,1),1:size(tmp,2),:,trialIx) = tmp;
        tmp =  aIM{trialIx};
        activIM(1:size(tmp,1),1:size(tmp,2),:,trialIx) = tmp;
    end
    clear mIM aIM tmp

    %Make template
    disp('Making template for aligning across trials...')
    maxshift = 5;
    M = squeeze(sum(meanIM, 3));
    samples = find(keepTrials(DMDix,:)); samples = samples(unique(round(linspace(1,length(samples),20))));
    template = makeConsensusTemplate(M(:,:,samples), maxshift);

    %align all mean images to template
    disp('Aligning across trials...')
    meanAligned = nan(size(meanIM),'single');
    actAligned = nan(size(meanIM,1),size(meanIM,2),1,nTrials,'single');
    corrCoeff = nan(1,nTrials);
    motOutput = nan(2,nTrials);
    Mpad = nan([size(template) size(M,3)],'single');
    Mpad(maxshift+(1:size(M,1)), maxshift+(1:size(M,2)),:) = M;

    fillval = min(template(:),[], 'omitnan')-1;
    tFFT = fft2(max(template, fillval));
    [rr,cc] = ndgrid(1:size(meanIM,1), 1:size(meanIM,2));
    for trialIx = nTrials:-1:1
        if ~keepTrials(DMDix,trialIx) || all(isnan(activIM(:,:,1,trialIx)), 'all')
            disp(['skipping trial, dmd:' int2str(trialIx) ' ' int2str(DMDix)])
            continue
        end
        disp(['trial: ' int2str(trialIx)])
        
        output1 = dftregistration_clipped(tFFT, fft2(max(Mpad(:,:,trialIx), fillval)),1,80);
        mot1 = [-output1(3) -output1(4)];
        [motOutput(:,trialIx), corrCoeff(trialIx)] = xcorr2_nans(Mpad(:,:,trialIx), template, round(mot1'), maxshift);

        for chIx = 1:size(meanIM,3)
            meanAligned(:,:,chIx,trialIx) = interp2(meanIM(:,:,chIx,trialIx), cc+motOutput(2,trialIx), rr+motOutput(1,trialIx));
        end
        actAligned(:,:,1,trialIx) = interp2(activIM(:,:,1,trialIx), cc+motOutput(2,trialIx), rr+motOutput(1,trialIx));
    end

    %identify outliers in alignment quality to determine valid trials
    corrThresh = min(0.90, median(corrCoeff, 'omitnan')-2*std(corrCoeff, 'omitmissing'));
    actValidPix = squeeze(mean(~isnan(actAligned(:,:,1,:)), [1 2]));
    validTrials= find(corrCoeff(:)>corrThresh & actValidPix(:)>mean(actValidPix)/2);
    if isempty(validTrials)
        error('SILo:NoValidAlignedTrials', ...
            ['No trials survived SILo cross-trial alignment QC for imaging path %d. ' ...
             'Inspect registration quality, discarded frames, and trial-table validity.'],DMDix);
    end
    exptSummary.meanIM{DMDix} = mean(meanAligned(:,:,:,validTrials),4, 'omitnan');

    %select source seeds on the aligned activity image using configured backend
    rawActIM = mean(actAligned(:,:,:,validTrials), 4, 'includenan');
    nanFrac = mean(isnan(actAligned(:,:,:,validTrials)), 4);
    % Mask out somata before either backend performs peak selection.
    somaMask = false(size(rawActIM));
    if ~isempty(ROIs)
        for rix = 1:numel(ROIs(DMDix).roiData)
            if contains(upper(ROIs(DMDix).roiData{rix}.Label), 'SOMA')
                tmp = ROIs(DMDix).roiData{rix}.mask;
                somaMask(1:size(tmp,1), 1:size(tmp,2)) = ...
                    somaMask(1:size(tmp,1), 1:size(tmp,2)) | tmp;
            end
        end
    end
    [sources, actIM, detectionInfo] = ...
        detectSourcesFromActivityImage(rawActIM, nanFrac, somaMask, params);
    exptSummary.actIM{DMDix} = actIM;
    exptSummary.rawActIM{DMDix} = rawActIM;
    exptSummary.sourceDetectionInfo{DMDix} = detectionInfo;
    sz = size(actIM);
    k = numel(sources.R);
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

    % Restore SILo-level parameters not owned by setParamsExtractTrial.
    params.sourceDetectionMethod = resolvedSourceDetectionMethod;
    params.savePerTrialSummary = resolvedSavePerTrialSummary;
    if ~isempty(resolvedMaxSynapseDensity)
        params.maxSynapseDensity = resolvedMaxSynapseDensity;
    end
    
    if isempty(ROIs) || isempty(ROIs(DMDix))
        roiData =[];
    else
        roiData = ROIs(DMDix).roiData;
    end

    if any(keepSources)
        highResTimer = tic;
        fns = trialTable.source_extraction.fn_raw(DMDix,:);
        if params.isSLAP2
            if isfield(trialTable, 'datadr') && ~isempty(trialTable.datadr)
                trialDataDr = trialTable.datadr;
            else
                trialDataDr = dr;
            end
            fls = trialTable.slap2_info.first_line(DMDix,:);
            els = trialTable.slap2_info.last_line(DMDix,:);
            E = processAllTrials_Async(trialDataDr, fns, fls, els, selPix, sources, discardFrames, alignData, meanAligned, motOutput, roiData, validTrials, params);
        else % non-SLAP2 (registered movies live under motion_correction)
            fls = cell(1,numel(fns));
            els = cell(1,numel(fns));
            E = processAllTrials_Async(mocodr, fns, fls, els, selPix, sources, discardFrames, alignData, meanAligned, motOutput, roiData, validTrials, params);
        end
        exptSummary.E(:,DMDix) = E;
        fprintf('Path %d high-resolution extraction: %.1f min\n', ...
            DMDix,toc(highResTimer)/60);
    end
    exptSummary.selPix{DMDix} = any(selPix,3);
    exptSummary.sources{DMDix} = sources;
    exptSummary.aData(:,DMDix) = alignData(:);
    exptSummary.userROIs{DMDix} = roiData;
    exptSummary.peaks{DMDix}= peaks;
    % Per-trial aligned images are only needed by per_trial_summary.h5.
    % When that optional output is disabled, release them immediately.
    exptSummary.perTrialMeanIMs{DMDix} = [];
    exptSummary.perTrialActIms{DMDix} = [];
    if params.savePerTrialSummary
        exptSummary.perTrialMeanIMsAligned{DMDix} = meanAligned;
        exptSummary.perTrialActIMsAligned{DMDix} = actAligned;
    else
        exptSummary.perTrialMeanIMsAligned{DMDix} = [];
        exptSummary.perTrialActIMsAligned{DMDix} = [];
    end
    exptSummary.perTrialAlignmentOffsets{DMDix} = motOutput;

    clear meanAligned meanIM activIM actAligned E
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
if params.savePerTrialSummary
    trySave(@() savePerTrialSummaryH5(fullfile(savedr, 'per_trial_summary.h5'), exptSummary, trialTable), 'per_trial_summary');
% Save detector-level diagnostics for direct SILo-vs-LoCo comparison.
benchmarkDiagnostics = struct();
benchmarkDiagnostics.sourceDetectionMethod = params.sourceDetectionMethod;
benchmarkDiagnostics.params = params;
benchmarkDiagnostics.rawActIM = exptSummary.rawActIM;
benchmarkDiagnostics.detectionIM = exptSummary.actIM;
benchmarkDiagnostics.sources = exptSummary.sources;
benchmarkDiagnostics.sourceDetectionInfo = exptSummary.sourceDetectionInfo;
diagnosticsName = ['source_detection_diagnostics_' params.sourceDetectionMethod '.mat'];
save(fullfile(savedr, diagnosticsName), 'benchmarkDiagnostics', '-v7.3');
else
    fprintf(['Skipping per_trial_summary.h5 (savePerTrialSummary=false). ' ...
        'Any existing file is left untouched.\n']);
end

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


fprintf('Total SILo runtime: %.2f hr\n',toc(runTimer)/3600);
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
