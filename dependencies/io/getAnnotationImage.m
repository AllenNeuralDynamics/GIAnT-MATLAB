function [IM, fn] = getAnnotationImage(trialTable, keepTrials, mocodr, DMDix)
%GETANNOTATIONIMAGE Mean image of the longest valid trial, the space ROIs are annotated in.
%
% Returns IM (H x W) and the registered file stem fn, or empty values if the
% imaging path has no valid trials. The longest trial (largest registered
% file) gives the cleanest mean image; short trials, often the first one, can
% be too noisy to annotate on or to register post-acquisition annotations to.
% Prefers meanIM from the alignment H5 and falls back to averaging the
% registered tiff.

IM = [];
fn = '';
validTrials = find(keepTrials(DMDix,:));
if isempty(validTrials)
    return
end
nBytes = zeros(size(validTrials));
for ix = 1:numel(validTrials)
    d = dir(fullfile(mocodr, trialTable.motion_correction.fn_reg_ds{DMDix, validTrials(ix)}));
    if ~isempty(d)
        nBytes(ix) = d(1).bytes;
    end
end
[~, best] = max(nBytes);
annTrial = validTrials(best);

[~, fn, ext] = fileparts(trialTable.motion_correction.fn_reg_ds{DMDix,annTrial});
adataFn = trialTable.motion_correction.fn_adata{DMDix, annTrial};
aDataAnnot = loadStructFromH5(fullfile(mocodr, adataFn));
if isfield(aDataAnnot, 'meanIM') && ~isempty(aDataAnnot.meanIM)
    IM = squeeze(mean(aDataAnnot.meanIM, 1, 'omitnan')); % numChannels x H x W -> H x W
else
    [IMtif, ~] = ScanImageTiffWrapper(fullfile(mocodr, [fn ext]));
    IM = squeeze(mean(IMtif, [3 4], 'omitnan'));
end
end
