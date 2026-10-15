function annotateROIs(dr_or_pathToTrialTable)

if ~nargin
    [trialTablefn, dr] =  uigetfile('*.h5', 'Select a trial_table file', '*trial_table*.h5' );
else
    %parse dr
    %_or_pathToTrialTable
    if exist(dr_or_pathToTrialTable, 'dir')
        dr = dr_or_pathToTrialTable;
        trialTablefn = 'trial_table.h5';
    else
        [dr trialTablefn ext] = fileparts(dr_or_pathToTrialTable);
        trialTablefn = [trialTablefn ext];
    end
end

%confirm that all files exist (also populates source_extraction.fn_raw)
[trialTable, keepTrials] = verifyFiles(trialTablefn, dr);
mocodr = fullfile(trialTable.savedr, 'motion_correction');
% for dmdKey = fieldnames(trialTable.slap2_info.ref_stack)'
%     trialTable.slap2_info.ref_stack.(dmdKey{1}).IM = []; %this uses a lot of memory and we won't need it
% end
nDMDs = size(trialTable.filename,1); % #imaging paths x #trials; non-SLAP2 is treated as 1 path
nTrials = size(trialTable.filename,2);

disp(['## ANNOTATING' newline 'Folder:'])
disp(dr)

%call up a GUI for the user to define Soma ROI and regions to exclude
%(uses existing annotations.h5, else post-acquisition annotations, else draws)
loadOrCreateAnnotations(trialTable, keepTrials, mocodr);
