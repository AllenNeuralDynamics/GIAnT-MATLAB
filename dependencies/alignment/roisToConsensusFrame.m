function roiEntry = roisToConsensusFrame(roiEntry, fnRegDS, motOutput, targetSz, warp)
%ROISTOCONSENSUSFRAME Move user ROIs from their annotation trial into SILo's consensus frame.
%
% roiEntry = roisToConsensusFrame(roiEntry, fnRegDS, motOutput, targetSz)
% roiEntry = roisToConsensusFrame(roiEntry, fnRegDS, motOutput, targetSz, warp)
%
% annotations.h5 stores ROIs in the frame of the registered trial they were
% drawn on (or registered to), roiEntry.fn. SILo aligns each trial's mean
% image to a consensus template, aligned(r,c) = trial(r + motOutput(1), c +
% motOutput(2)), and uses that frame for soma exclusion and ROI extraction.
% The ROIs are therefore shifted by -motOutput(:, annotation trial).
%
% With warp (SILo's across-trial alignment, see trialWarpCoords), an affine or
% patchwise map of the annotation trial is applied in full instead: a mask
% pixel x of the consensus frame is set where the trial's mask is set at
% trialWarpCoords(warp, trial, x) (nearest pixel), and ROI geometry is mapped
% with trialWarpInverse. A shift of the annotation trial is applied as above.
%
% fnRegDS   - registered file names for this imaging path (1 x nTrials)
% motOutput - SILo's cross-trial offsets for this path (2 x nTrials)
% targetSz  - [rows cols] of the consensus frame
% warp      - optional; SILo's across-trial maps for this path
%
% If the annotation trial cannot be identified, or was not aligned (NaN
% offset), the ROIs are not shifted, with a warning. Their masks are still
% cropped or padded to targetSz: soma exclusion and the experiment summary
% index them on the consensus frame and fail on any other size.

if ~isfield(roiEntry, 'roiData') || isempty(roiEntry.roiData)
    return
end
if nargin < 5
    warp = [];
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
elseif ~isempty(warp) && ~isShift(warp, trialIx)
    roiEntry.roiData = warpRoiData(roiEntry.roiData, warp, trialIx, targetSz);
    return
else
    offset = motOutput(:, trialIx);
end
roiEntry.roiData = shiftRoiData(roiEntry.roiData, -offset(:)', targetSz);
end


function tf = isShift(warp, t)
tf = isempty(warp.nodes) && isequal(warp.affine(:,1:2,t), eye(2));
end


function roiData = warpRoiData(roiData, warp, t, targetSz)
%masks: each consensus-frame pixel takes the trial mask's value at its mapped
%position; geometry ([x y], images.roi order) is mapped from the trial's image
[cg, rg] = meshgrid(1:targetSz(2), 1:targetSz(1));
[rT, cT] = trialWarpCoords(warp, t, rg, cg);
rT = round(rT); cT = round(cT);
for rix = 1:numel(roiData)
    S = roiData{rix};
    if isfield(S, 'mask') && ~isempty(S.mask)
        m = logical(S.mask);
        in = rT >= 1 & rT <= size(m, 1) & cT >= 1 & cT <= size(m, 2);
        out = false(targetSz);
        out(in) = m(sub2ind(size(m), rT(in), cT(in)));
        S.mask = out;
    end
    for fnm = {'Center', 'Position', 'Vertices'}
        if isfield(S, fnm{1}) && ~isempty(S.(fnm{1}))
            xy = double(S.(fnm{1}));
            [r, c] = trialWarpInverse(warp, t, xy(:,2), xy(:,1));
            S.(fnm{1}) = [c(:) r(:)];
        end
    end
    roiData{rix} = S;
end
end


function s = stemOf(f)
[~, s] = fileparts(char(f));
end
