function trialTable = buildTrialTableSLAP2(dr,savedr,useAllFiles)
%BUILDTRIALTABLESLAP2 Organize multi-trial SLAP2 recordings and reference images.
%   Builds the trial_table struct documented in README.md and writes it to
%   fullfile(savedr, 'trial_table.h5'). SLAP2-specific metadata (reference
%   stacks, first/last lines, trial timing) is nested under slap2_info.
%
%   This function addresses a bug in SLAP2 trial numbering as of Feb 2024,
%   where trial numbers sometimes fail to increment. That makes some files
%   extra long, and subsequent trial numbers get out of sync.
%
%   Each path's files are put in acquisition order by the timestamp and
%   trial number in their names, not by modification time, which copies
%   (e.g. from S3) don't preserve. Trial times come from the FPGA line
%   timestamps when the reader and .meta file provide them, and otherwise
%   from file modification times.

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
    numLines = zeros(1, numel(files));
    firstLineTime = nan(1, numel(files));
    for fIx = length(files):-1:1
        hDat = slap2.Slap2DataFile([dr filesep files(fIx).name]);
        numLines(fIx) = hDat.totalNumLines;
        firstLineTime(fIx) = readFirstLineTime(hDat);
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

    % Put each path's files in acquisition order. File modification times
    % can't be trusted for this: copies (e.g. downloads from S3) carry the
    % copy time, so use the acquisition timestamp and trial number in the
    % file names instead.
    DMD1files = files(isDMD1);
    DMD2files = files(~isDMD1);
    [order1, namesOk1] = acquisitionOrder({DMD1files.name});
    [order2, namesOk2] = acquisitionOrder({DMD2files.name});
    namesOk = namesOk1 && namesOk2;
    if ~namesOk
        warning('buildTrialTableSLAP2:FileOrderFromModTime', ...
            ['Could not read acquisition times and trial numbers from the .dat file names ' ...
            'in epoch %d; ordering files by modification time instead, which is wrong ' ...
            'for copied files.'], eIx);
        [~, order1] = sort([DMD1files.datenum], 'ascend');
        [~, order2] = sort([DMD2files.datenum], 'ascend');
    end
    DMD1files = DMD1files(order1);
    numLinesDMD1 = numLines(isDMD1); numLinesDMD1 = numLinesDMD1(order1);
    startTimeDMD1 = firstLineTime(isDMD1); startTimeDMD1 = startTimeDMD1(order1);

    DMD2files = DMD2files(order2);
    numLinesDMD2 = numLines(~isDMD1); numLinesDMD2 = numLinesDMD2(order2);
    startTimeDMD2 = firstLineTime(~isDMD1); startTimeDMD2 = startTimeDMD2(order2);

    % Trial times come from the FPGA line clock when the reader provides it.
    % Otherwise fall back to file modification times, and warn if those look
    % like copy times rather than the times the files were written.
    useLineClock = all(isfinite(firstLineTime));
    if ~useLineClock
        reason = copyTimeReason(DMD1files, numLinesDMD1, linePeriod_s, namesOk);
        if isempty(reason)
            reason = copyTimeReason(DMD2files, numLinesDMD2, linePeriod_s, namesOk);
        end
        if ~isempty(reason)
            warning('buildTrialTableSLAP2:ModTimesNotAcquisitionTimes', ...
                ['File modification times in epoch %d do not look like acquisition times (%s). ' ...
                'trial_end_time_from_pc and trial_start_time_inferred come from these times ' ...
                'and are probably copy times.'], eIx, reason);
        end
    end
    lineDays = linePeriod_s/86400; % line period in datenum units

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

                if useLineClock
                    slap2_info.trial_end_time_from_pc(trueTrialIx) = startTimeDMD1(fileIx) + (trialEdges(trialIx+1)-1)*lineDays;
                    slap2_info.trial_start_time_inferred(trueTrialIx) = startTimeDMD1(fileIx) + (trialEdges(trialIx)-1)*lineDays;
                else
                    slap2_info.trial_end_time_from_pc(trueTrialIx) = DMD1files(fileIx).datenum - datenum(duration(0,0,(numLinesDMD1(fileIx)-trialEdges(trialIx+1)).*linePeriod_s));
                    slap2_info.trial_start_time_inferred(trueTrialIx) = DMD1files(fileIx).datenum - datenum(duration(0,0,(numLinesDMD1(fileIx)-trialEdges(trialIx)).*linePeriod_s));
                end
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

        % At least one path starts a new file at every trial, and that
        % file's first line timestamp is the trial's start.
        if accumLines(1) == 0
            lineClockStart = startTimeDMD1(lastDMD1fIx+1);
        else
            lineClockStart = startTimeDMD2(lastDMD2fIx+1);
        end

        if abs(cumLines1(1)-cumLines2(1))<lineDiffThresh
            endTimeFromMtime = DMD1files(lastDMD1fIx+1).datenum;

            lastDMD1fIx = lastDMD1fIx+1; %we finished this DMD1 file
            lastDMD2fIx = lastDMD2fIx+1; %we finished this DMD2 file
            accumLines = [0 0]; %reset accumulated lines

        elseif cumLines1(1) < cumLines2(1)
            endTimeFromMtime = DMD1files(lastDMD1fIx+1).datenum;

            lastDMD1fIx = lastDMD1fIx+1; %we finished this DMD1 file
            accumLines(1) = 0;
            accumLines(2) = accumLines(2) + nLines;
        elseif cumLines2(1) < cumLines1(1)
            endTimeFromMtime = DMD2files(lastDMD2fIx+1).datenum;

            lastDMD2fIx = lastDMD2fIx+1; %we finished this DMD2 file
            accumLines(1) = accumLines(1) + nLines;
            accumLines(2) = 0;
        end

        if useLineClock
            slap2_info.trial_start_time_inferred(trueTrialIx) = lineClockStart;
            slap2_info.trial_end_time_from_pc(trueTrialIx) = lineClockStart + nLines*lineDays;
        else
            slap2_info.trial_end_time_from_pc(trueTrialIx) = endTimeFromMtime;
            slap2_info.trial_start_time_inferred(trueTrialIx) = endTimeFromMtime - datenum(duration(0,0,nLines.*linePeriod_s));
        end
    end
    end
end

trialTable.slap2_info = slap2_info;

saveStructToH5(trialTable, fullfile(savedr, 'trial_table.h5'));
end

function [order, ok] = acquisitionOrder(names)
%ACQUISITIONORDER Order SLAP2 .dat files by acquisition time, then trial.
%   Names look like activity_20250314_111223_DMD1-TRIAL000002.dat; files
%   from continuous acquisitions end in -CYCLE-000000 instead of, or after,
%   the trial number. ok is false when files from different acquisitions
%   have no timestamp in their names to tell them apart, or when two files
%   get the same position (e.g. no trial numbers). order is then the
%   identity.
n = numel(names);
order = 1:n;
acqTime = nan(n, 1);
trial = zeros(n, 1);
base = cell(n, 1);
for k = 1:n
    stem = regexprep(names{k}, '-CYCLE-?\d+\.dat$|\.dat$', '', 'ignorecase');
    trialTok = regexpi(stem, '-TRIAL(\d+)$', 'tokens', 'once');
    if ~isempty(trialTok)
        trial(k) = str2double(trialTok{1});
    end
    base{k} = regexprep(stem, '-TRIAL\d+$', '', 'ignorecase');
    stamp = regexp(base{k}, '(?<!\d)(\d{8})_(\d{6})(?!\d)', 'tokens');
    if ~isempty(stamp)
        acqTime(k) = parseNameTimestamp([stamp{end}{:}]);
    end
end
if all(isfinite(acqTime))
    key = [acqTime trial];
elseif numel(unique(base)) == 1
    key = trial; % a single acquisition: the trial number alone sets the order
else
    ok = false;
    return
end
ok = size(unique(key, 'rows'), 1) == n;
if ok
    [~, order] = sortrows(key);
    order = order';
end
end

function t = parseNameTimestamp(str)
% datenum of a yyyyMMddHHmmss string, or NaN if it isn't a valid date.
v = sscanf(str, '%4d%2d%2d%2d%2d%2d')';
t = NaN;
if numel(v) == 6 && v(2) >= 1 && v(2) <= 12 && v(3) >= 1 && v(3) <= 31 ...
        && v(4) < 24 && v(5) < 60 && v(6) < 60
    t = datenum(v);
end
end

function t = readFirstLineTime(hDat)
%READFIRSTLINETIME Acquisition time (datenum) of a .dat file's first line.
%   SLAP2 stamps every line with the FPGA clock, and the .meta file's
%   fpgaTimeReference is the time of FPGA tick 0. NaN if the reader or the
%   metadata don't provide these.
t = NaN;
if isprop(hDat, 'hMultiDataFiles')
    src = hDat.hMultiDataFiles;
elseif isprop(hDat, 'hDataFile')
    src = hDat.hDataFile;
else
    return
end
ticks = getMember(src, 'firstLineTimestamp');
clock_Hz = getMember(getMember(src, 'header'), 'fpgaSystemClock_Hz');
ref = getMember(getMember(src, 'metaData'), 'fpgaTimeReference');
if isscalar(ref) && isdatetime(ref) && ~isnat(ref) && isscalar(ticks) && ticks > 0 ...
        && isscalar(clock_Hz) && clock_Hz > 0
    t = datenum(ref + seconds(double(ticks)/double(clock_Hz)));
end
end

function v = getMember(x, name)
% Field of a struct or property of an object; [] if it has none.
v = [];
if (isstruct(x) && isscalar(x) && isfield(x, name)) || (isobject(x) && isprop(x, name))
    v = x.(name);
end
end

function reason = copyTimeReason(files, numLines, linePeriod_s, checkOrder)
%COPYTIMEREASON Why files' modification times look like copy times, or ''.
%   SLAP2 writes each .dat file as its trial ends, so in acquisition order
%   the modification times increase and span at least the recording time of
%   all files but one. Copies carry the copy time instead, often all within
%   a few seconds and in arbitrary order. checkOrder says whether files are
%   in acquisition order (from their names) rather than in mtime order.
reason = '';
if numel(files) < 2
    return
end
mtime_s = [files.datenum]*86400;
recorded_s = numLines*linePeriod_s;
minSpan_s = sum(recorded_s) - max(recorded_s);
span_s = max(mtime_s) - min(mtime_s);
if checkOrder && any(diff(mtime_s) < -2) % allow for coarse mtime resolution
    reason = 'they are out of order relative to the file names';
elseif span_s < 0.5*minSpan_s
    reason = sprintf('they span %.0f s, but the trials add up to at least %.0f s', ...
        span_s, minSpan_s);
end
end
