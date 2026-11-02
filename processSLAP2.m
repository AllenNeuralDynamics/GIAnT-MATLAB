function [trialTable, mocoParams, siloParams] = processSLAP2(dr, savedr, useAllFiles)
%PROCESSSLAP2 Run the GIAnT SLAP2 pipeline from a single call.
%
%   processSLAP2()
%       Select a SLAP2 root folder interactively. Outputs are written to
%       that folder.
%
%   processSLAP2(dr)
%       Process the SLAP2 dataset rooted at dr. For the standard AIND
%       layout, dr should be the parent "slap2" folder containing
%       dynamic_data/ and static_data/.
%
%   processSLAP2(dr, savedr)
%       Write GIAnT outputs to savedr instead of dr.
%
%   processSLAP2(dr, savedr, useAllFiles)
%       If useAllFiles is true, buildTrialTableSLAP2 places all discovered
%       acquisitions into one epoch. The default is false, which invokes
%       the normal GIAnT epoch-selection dialog.
%
%   [trialTable, mocoParams, siloParams] = processSLAP2(...)
%       Also returns the final trial table and the parameter structures used
%       for motion correction and source extraction.
%
% Pipeline:
%   1. buildTrialTableSLAP2
%   2. MultiRoiRegistration
%   3. SILo
%
% The GIAnT optionsGUI is shown once for MultiRoiRegistration and once for
% SILo. SILo is forced to run with isSLAP2 = true. If drawUserRois is true,
% SILo will also open its normal annotation GUI during source extraction.

% Put this whole repo ahead of anything else on the MATLAB path. Other
% checkouts (e.g. SLAP2-Voltage-Extraction) carry their own copies of
% MultiRoiRegistration, SILo and several dependencies; only the top-level
% folders were previously left to path luck. Parallel workers are pinned the
% same way by refreshPoolWorkerCode when a pool is created or reused.
here = fileparts(mfilename('fullpath'));
repoRoot = here;
while ~isfolder(fullfile(repoRoot,'dependencies')) && ~strcmp(repoRoot, fileparts(repoRoot))
    repoRoot = fileparts(repoRoot);   % walk up until a 'dependencies' folder is found
end
assert(isfolder(fullfile(repoRoot,'dependencies')), 'No dependencies folder found above %s', here);
addpath(genpath(repoRoot), '-begin');
rehash;

%% Resolve input/output directories
if nargin < 1 || isempty(dr)
    dr = uigetdir(pwd, 'Select SLAP2 root directory');
    if isequal(dr, 0)
        error('processSLAP2:Cancelled', 'SLAP2 directory selection cancelled.');
    end
end

if nargin < 2 || isempty(savedr)
    savedr = dr;
end

if nargin < 3 || isempty(useAllFiles)
    useAllFiles = false;
end

dr = char(dr);
savedr = char(savedr);
useAllFiles = logical(useAllFiles);

if ~isfolder(dr)
    error('processSLAP2:InvalidDataDirectory', ...
        'SLAP2 data directory does not exist:\n%s', dr);
end

if ~isfolder(savedr)
    mkdir(savedr);
end

%% Check required GIAnT functions before starting
requiredFunctions = {
    'buildTrialTableSLAP2'
    'MultiRoiRegistration'
    'SILo'
    'setParams'
    'loadStructFromH5'
    'extractTrial'
    'processAllTrials_Async'
    'getActImPeaks'
    'interpFramesSelectedBatch'
    'interpFrameTranslationChannels'
    'xcorr2_nans_weighted_fast'
    'Fast_BigTiff_Write'
    'slap2.Slap2DataFile'
    'fmincon'
    'lsqlin'
    'imgaussfilt'
    'medfilt1'
    'prctile'
};

for i = 1:numel(requiredFunctions)
    if isempty(which(requiredFunctions{i}))
        error('processSLAP2:MissingDependency', ...
            'Required function is not on the MATLAB path: %s', ...
            requiredFunctions{i});
    end
end

fprintf('\n============================================================\n');
fprintf('GIAnT SLAP2 processing\n');
fprintf('Data root:   %s\n', dr);
fprintf('Output root: %s\n', savedr);
fprintf('============================================================\n\n');

%% 1. Build or reuse trial table
trialTablePath = fullfile(savedr, 'trial_table.h5');

if isfile(trialTablePath)
    choice = questdlg( ...
        sprintf('A trial_table.h5 already exists in:\n%s\n\nWhat would you like to do?', savedr), ...
        'Existing GIAnT trial table', ...
        'Reuse', 'Rebuild', 'Cancel', 'Reuse');

    switch choice
        case 'Reuse'
            fprintf('Reusing existing trial table:\n  %s\n\n', trialTablePath);
            trialTable = loadStructFromH5(trialTablePath);

        case 'Rebuild'
            fprintf('Rebuilding trial table...\n');
            trialTable = buildTrialTableSLAP2(dr, savedr, useAllFiles);

        otherwise
            error('processSLAP2:Cancelled', 'Processing cancelled by user.');
    end
else
    fprintf('Building trial table...\n');
    trialTable = buildTrialTableSLAP2(dr, savedr, useAllFiles);
end

% Basic guard against an incomplete trial table.
requiredTrialFields = {'datadr', 'savedr', 'filename', 'true_trial_ix', 'epoch', 'slap2_info'};
for i = 1:numel(requiredTrialFields)
    if ~isfield(trialTable, requiredTrialFields{i})
        error('processSLAP2:IncompleteTrialTable', ...
            'Trial table is missing required field "%s".', requiredTrialFields{i});
    end
end

if isempty(trialTable.filename)
    error('processSLAP2:EmptyTrialTable', ...
        'The trial table contains no acquisition trials.');
end

if ~isfield(trialTable.slap2_info, 'first_line') || ...
        ~isfield(trialTable.slap2_info, 'last_line')
    error('processSLAP2:IncompleteTrialTable', ...
        'Trial table is missing SLAP2 first_line/last_line information.');
end

fprintf('Trial table ready: %d imaging path(s) x %d analysis trial(s)\n\n', ...
    size(trialTable.filename, 1), size(trialTable.filename, 2));

%% 2. Configure motion correction
% Calling setParams directly shows optionsGUI here. Passing the resulting
% struct to MultiRoiRegistration prevents the function from opening a
% second, duplicate optionsGUI.
fprintf('Configure MultiRoiRegistration options...\n');
mocoParams = setParams('MultiRoiRegistration');

%% 3. Configure source extraction
% Seed the GUI with isSLAP2=true, because setParams defaults this to false.
% forceGUI=true preserves the interactive optionsGUI while using the SLAP2
% default appropriate for this wrapper.
fprintf('Configure SILo source-extraction options...\n');
siloParams = setParams('SILo', struct('isSLAP2', true), true);

% processSLAP2 is specifically a SLAP2 wrapper, so do not permit this flag
% to be accidentally disabled in the GUI.
siloParams.isSLAP2 = true;

fprintf('\nConfiguration complete. Starting processing...\n\n');

%% 4. Motion correction
fprintf('============================================================\n');
fprintf('STEP 1/2: MultiRoiRegistration\n');
fprintf('============================================================\n');

mocoParams = MultiRoiRegistration(trialTablePath, mocoParams);

%% 5. Source extraction
fprintf('\n============================================================\n');
fprintf('STEP 2/2: SILo source extraction\n');
fprintf('============================================================\n');

siloParams = SILo(trialTablePath, siloParams);

%% Reload final trial table so returned value includes downstream metadata
trialTable = loadStructFromH5(trialTablePath);

fprintf('\n============================================================\n');
fprintf('GIAnT SLAP2 processing complete.\n');
fprintf('Trial table:       %s\n', trialTablePath);
fprintf('Motion correction: %s\n', fullfile(savedr, 'motion_correction'));
fprintf('Source extraction: %s\n', fullfile(savedr, 'source_extraction'));
fprintf('============================================================\n\n');

end
