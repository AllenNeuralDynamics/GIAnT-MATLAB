function roiEntry = roisToConsensusFrame(roiEntry, fnRegDS, motOutput, targetSz)
%ROISTOCONSENSUSFRAME Move user ROIs from their annotation trial into SILo's consensus frame.
%
% roiEntry = roisToConsensusFrame(roiEntry, fnRegDS, motOutput, targetSz)
%
% annotations.h5 stores ROIs in the frame of the registered trial they were
% drawn on (or registered to), roiEntry.fn. SILo aligns each trial's mean
% image to a consensus template, aligned(r,c) = trial(r + motOutput(1), c +
% motOutput(2)), and uses that frame for soma exclusion and ROI extraction.
% The ROIs are therefore shifted by -motOutput(:, annotation trial).
%
% fnRegDS   - registered file names for this imaging path (1 x nTrials)
% motOutput - SILo's cross-trial offsets for this path (2 x nTrials)
% targetSz  - [rows cols] of the consensus frame
%
% If the annotation trial cannot be identified, or was not aligned (NaN
% offset), the ROIs are not shifted, with a warning. Their masks are still
% cropped or padded to targetSz: soma exclusion and the experiment summary
% index them on the consensus frame and fail on any other size.

if ~isfield(roiEntry, 'roiData') || isempty(roiEntry.roiData)
    return
end
%fn is normally the registered file's stem (getAnnotationImage), which can itself
%contain dots, so match it whole before trying it without an extension
[~, name, ext] = fileparts(char(roiEntry.fn));
annName = [name ext];
stems = cellfun(@(f) stemOf(f), fnRegDS, 'UniformOutput', false);
trialIx = find(strcmp(stems, annName), 1);
if isempty(trialIx)
    trialIx = find(strcmp(stems, name), 1);
end
offset = [0; 0];
if isempty(annName) || isempty(trialIx)
    warning('roisToConsensusFrame:UnknownTrial', ...
        'Annotation trial "%s" is not in the trial table; using its ROIs without cross-trial alignment.', annName);
elseif any(~isfinite(motOutput(:, trialIx)))
    warning('roisToConsensusFrame:TrialNotAligned', ...
        'Annotation trial "%s" was not aligned across trials; using its ROIs without cross-trial alignment.', annName);
else
    offset = motOutput(:, trialIx);
end
roiEntry.roiData = shiftRoiData(roiEntry.roiData, -offset(:)', targetSz);
end


function s = stemOf(f)
[~, s] = fileparts(char(f));
end
