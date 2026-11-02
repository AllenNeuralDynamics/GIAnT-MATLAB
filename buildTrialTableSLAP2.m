function trialTable = buildTrialTableSLAP2(dr,savedr,useAllFiles)
%BUILDTRIALTABLESLAP2 Organize multi-trial SLAP2 recordings and reference images.
%   Builds the trial_table struct documented in README.md and writes it to
%   fullfile(savedr, 'trial_table.h5'). SLAP2-specific metadata (reference
%   stacks, first/last lines, trial timing) is nested under slap2_info.
%
%   This function addresses a bug in SLAP2 trial numbering as of Feb 2024,
%   where trial numbers sometimes fail to increment. That makes some files
%   extra long, and subsequent trial numbers get out of sync.

%parameters
lineDiffThresh = 2000; %difference threshold for calling two recordings the same length, in lines. ~0.2 seconds
multiCycleLinesPerTrial = 200000; % Break up continuous acquisitoins into blocks of this many lines

%get a list of dat files in a given folder
if ~nargin
    dr = uigetdir;
end
if nargin < 2
    savedr = dr;
end
if nargin < 3
    useAllFiles = false;
end
unpickedfiles = dir([dr filesep '*.dat']);

%remove extra 'multicycle' files from list; they will be represented by first file
pattern = 'CYCLE-?(\d+)'; 
removeFiles = false(1,numel(unpickedfiles));
for ix = 1:length(unpickedfiles)
        % Use regexp to find matches and extract the number part
        tokens = regexp(unpickedfiles(ix).name, pattern, 'tokens');      
        if ~isempty(tokens) && str2double(tokens{1}{1})>0
            removeFiles(ix) = true;
        end
end
unpickedfiles = unpickedfiles(~removeFiles);

if useAllFiles
    epoch = 1;
    epochfiles{1} = unpickedfiles;
else
    epoch = 0;
    while ~isempty(unpickedfiles)
        [indx,tf] = listdlg('ListString',{unpickedfiles.name}, 'PromptString',['Select files for EPOCH ' int2str(epoch)]);
        if ~tf
            break
        end
        epoch = epoch+1;
        epochfiles{epoch} = unpickedfiles(indx);
        unpickedfiles(indx) = [];
    end
end

%first, load the reference images, deduce the imaging plane of the ROI, and
% %compute the soma ROI
DMDixs = [1 2];
nDMDs = numel(DMDixs);
refStack = struct();
%find the files within the entire folder structure named REFERENCE
for DMDix = DMDixs
    pathKey = ['Path' int2str(DMDix)];
    list = dir([dr filesep '**' filesep '*DMD' int2str(DMDix) '*_CONFIG1-REFERENCE*']);
    if isempty(list)
        list = dir([dr filesep '**' filesep '*DMD' int2str(DMDix) '*-REFERENCE*']);
    end
    switch length(list)
    case 0
        error(['Could not find reference image within folder: ' dr])
    case 1
        DMD1refFn = fullfile(list(1).folder, list(1).name);
        [refStack.(pathKey).IM, ~, desc] = ScanImageTiffWrapper(DMD1refFn);

        %load metadata for the reference stack and BCI ROI
        zIx = 0; Zs = []; channels = [];
        for imIx = 1:numel(desc)
                jj = jsondecode(desc{imIx});
                if ~any(channels == jj.channel)
                    channels = [channels jj.channel];
                end
                if imIx==1
                    accumChan = jj.channel;
                end
                if jj.channel == accumChan
                    zIx = zIx +1;
                    Zs(zIx) = jj.z;
                end
        end
        refStack.(pathKey).channels = channels;
        refStack.(pathKey).Zs = Zs;
        refStack.(pathKey).dmdPixel2SampleTransform = jj.dmdPixel2SampleTransform;

    otherwise
        list = dir([dr filesep '**' filesep '*DMD' int2str(DMDix) '_CONFIG1-REFERENCE*.tif']);
        switch length(list)
            case 0
                error('Too many reference stacks found in the specified directory!');
            case 1
                DMD1refFn = fullfile(list(1).folder, list(1).name);
                [refStack.(pathKey).IM, ~, desc] = ScanImageTiffWrapper(DMD1refFn);

                %load metadata for the reference stack and BCI ROI
                zIx = 0; Zs = []; channels = [];
                for imIx = 1:numel(desc)
                        jj = jsondecode(desc{imIx});
                        if ~any(channels == jj.channel)
                            channels = [channels jj.channel];
                        end
                        if imIx==1
                            accumChan = jj.channel;
                        end
                        if jj.channel == accumChan
                            zIx = zIx +1;
                            Zs(zIx) = jj.z;
                        end
                end
                refStack.(pathKey).channels = channels;
                refStack.(pathKey).Zs = Zs;
                refStack.(pathKey).dmdPixel2SampleTransform = jj.dmdPixel2SampleTransform;
            otherwise
                error('Too many CONFIG1 reference stacks found in the specified directory!');
        end
    end


end

trialTable.datadr = dr;
trialTable.savedr = savedr;

trialTable.filename = {};

slap2_info = struct();
slap2_info.ref_stack = refStack; clear refStack;
slap2_info.first_line = [];
slap2_info.last_line = [];
slap2_info.trial_end_time_from_pc = [];
slap2_info.trial_start_time_inferred = [];

trueTrialIx = 0;
for eIx = 1:epoch %for each epoch
    %get number of lines from dat file timestamp
    files = epochfiles{eIx};
    disp(['Loading metadata from ' int2str(length(epochfiles{eIx})) ' DAT files...']);
    for fIx = length(files):-1:1
        hDat = slap2.Slap2DataFile([dr filesep files(fIx).name]);
        numLines(fIx) = hDat.totalNumLines;
    end

    if isprop(hDat, 'hDataFile')
        linePeriod_s = hDat.hDataFile.metaData.linePeriod_s;
    elseif isprop(hDat, 'hMultiDataFiles')
        linePeriod_s = hDat.hMultiDataFiles.metaData.linePeriod_s;
    else
        error('Could not read data file metaData')
    end

    disp('done loading ')
    isDMD1 = cellfun(@(x)(~isempty(x)), strfind({files.name}, 'DMD1'));

    DMD1files = files(isDMD1);
    [~, sortorder] = sort([DMD1files.datenum], 'ascend');
    DMD1files = DMD1files(sortorder);
    numLinesDMD1 = numLines(isDMD1); numLinesDMD1 = numLinesDMD1(sortorder);

    DMD2files = files(~isDMD1);
    [~, sortorder] = sort([DMD2files.datenum], 'ascend');
    DMD2files = DMD2files(sortorder);
    numLinesDMD2 = numLines(~isDMD1); numLinesDMD2 = numLinesDMD2(sortorder);

    %Each epoch should consist of only continuous or only triggered
    %acquisitions
    if contains(DMD1files(1).name, '-CYCLE') %Continuous acquisition mode
        assert(numel(DMD1files)==numel(DMD2files), 'There was an unequal amount of DMD1 and DMD2 files for an epoch using continuous acquisitions');
        for fileIx = 1:numel(DMD1files)
            nLinesTot = min(numLinesDMD1(fileIx), numLinesDMD2(fileIx));
            nTrialsInFile = ceil(nLinesTot/multiCycleLinesPerTrial);
            trialEdges = linspace(1,nLinesTot+1, nTrialsInFile+1);

            for trialIx = 1:nTrialsInFile
                trueTrialIx = trueTrialIx+1;
                trialTable.filename{1,trueTrialIx} = DMD1files(fileIx).name;
                trialTable.filename{2,trueTrialIx} = DMD2files(fileIx).name;
                slap2_info.first_line(1,trueTrialIx) = trialEdges(trialIx);
                slap2_info.first_line(2,trueTrialIx) = trialEdges(trialIx);
                slap2_info.last_line(1,trueTrialIx) = trialEdges(trialIx+1);
                slap2_info.last_line(2,trueTrialIx) = trialEdges(trialIx+1);
                trialTable.true_trial_ix(1:nDMDs, trueTrialIx) = trueTrialIx;
                trialTable.epoch(1:nDMDs, trueTrialIx) = eIx;

                slap2_info.trial_end_time_from_pc(trueTrialIx) = DMD1files(fileIx).datenum - datenum(duration(0,0,(numLinesDMD1(fileIx)-trialEdges(trialIx+1)).*linePeriod_s));
                slap2_info.trial_start_time_inferred(trueTrialIx) = DMD1files(fileIx).datenum - datenum(duration(0,0,(numLinesDMD1(fileIx)-trialEdges(trialIx)).*linePeriod_s));
            end
        end
    else %triggered acquisition more+
    lastDMD1fIx = 0;
    lastDMD2fIx = 0;
    accumLines = [0 0];
    while lastDMD2fIx<length(DMD2files) && lastDMD1fIx<length(DMD1files)

        %confirm that the trials match up AT SOME POINT SOON
        cumLines1 = cumsum(numLinesDMD1(lastDMD1fIx+1: min(end, lastDMD1fIx+5)))-accumLines(1);
        cumLines2 = cumsum(numLinesDMD2(lastDMD2fIx+1: min(end, lastDMD2fIx+5)))-accumLines(2);
        assert(min(abs(cumLines1 - cumLines2'), [],'all')<lineDiffThresh, 'Error lining up trials across DMDs!')

        nLines = min(cumLines1(1), cumLines2(1));
        trueTrialIx = trueTrialIx+1;
        trialTable.filename{1,trueTrialIx} = DMD1files(lastDMD1fIx+1).name;
        trialTable.filename{2,trueTrialIx} = DMD2files(lastDMD2fIx+1).name;
        slap2_info.first_line(1,trueTrialIx) = accumLines(1)+1;
        slap2_info.first_line(2,trueTrialIx) = accumLines(2)+1;
        slap2_info.last_line(1,trueTrialIx) = accumLines(1)+nLines;
        slap2_info.last_line(2,trueTrialIx) = accumLines(2)+nLines;
        trialTable.true_trial_ix(1:nDMDs, trueTrialIx) = trueTrialIx;
        trialTable.epoch(1:nDMDs, trueTrialIx) = eIx;

        if abs(cumLines1(1)-cumLines2(1))<lineDiffThresh
            slap2_info.trial_end_time_from_pc(trueTrialIx) = DMD1files(lastDMD1fIx+1).datenum;
            slap2_info.trial_start_time_inferred(trueTrialIx) = slap2_info.trial_end_time_from_pc(trueTrialIx) - datenum(duration(0,0,nLines.*linePeriod_s));

            lastDMD1fIx = lastDMD1fIx+1; %we finished this DMD1 file
            lastDMD2fIx = lastDMD2fIx+1; %we finished this DMD2 file
            accumLines = [0 0]; %reset accumulated lines

        elseif cumLines1(1) < cumLines2(1)
            slap2_info.trial_end_time_from_pc(trueTrialIx) = DMD1files(lastDMD1fIx+1).datenum;
            slap2_info.trial_start_time_inferred(trueTrialIx) = slap2_info.trial_end_time_from_pc(trueTrialIx) - datenum(duration(0,0,nLines.*linePeriod_s));

            lastDMD1fIx = lastDMD1fIx+1; %we finished this DMD1 file
            accumLines(1) = 0;
            accumLines(2) = accumLines(2) + nLines;
        elseif cumLines2(1) < cumLines1(1)
            slap2_info.trial_end_time_from_pc(trueTrialIx) = DMD2files(lastDMD2fIx+1).datenum;
            slap2_info.trial_start_time_inferred(trueTrialIx) = slap2_info.trial_end_time_from_pc(trueTrialIx) - datenum(duration(0,0,nLines.*linePeriod_s));

            lastDMD2fIx = lastDMD2fIx+1; %we finished this DMD1 file
            accumLines(1) = accumLines(1) + nLines;
            accumLines(2) = 0;
        end
    end
    end
end

trialTable.slap2_info = slap2_info;

saveStructToH5(trialTable, fullfile(savedr, 'trial_table.h5'));
end
