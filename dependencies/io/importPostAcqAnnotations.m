function [ROIs, imported] = importPostAcqAnnotations(trialTable, keepTrials, mocodr)
%IMPORTPOSTACQANNOTATIONS Align post-acquisition annotations to the registered data.
%
% [ROIs, imported] = importPostAcqAnnotations(trialTable, keepTrials, mocodr)
%
% Looks for postacq_annotations.h5 (written by postAcquisitionAnnotation) in
% trialTable.datadr or its subfolders. ROIs in that file were drawn on raw
% acquisition images, so they are translated into the space of each imaging
% path's annotation image (see getAnnotationImage) before use.
% imported(DMDix) is true for paths whose annotations were registered; other
% paths (not annotated post-acquisition, or registration failed) have empty
% roiData. imported is all false if no post-acquisition file is found.

nDMDs = size(trialTable.filename, 1);
imported = false(1, nDMDs);
ROIs = repmat(struct('dr', mocodr, 'fn', '', 'roiData', {{}}, ...
    'alignmentShiftRC', [], 'alignmentNCC', []), 1, nDMDs);
fnPostAcq = findPostAcqAnnotationFile(trialTable.datadr);
if isempty(fnPostAcq)
    return
end
disp(['Importing post-acquisition annotations: ' fnPostAcq])
[postAcq, schemaVersion] = loadAnnotationsH5(fnPostAcq);
if schemaVersion ~= annotationsSchemaVersion()
    %ROI geometry and the image orientation they were drawn in must match
    %what this version expects, or ROIs would be placed incorrectly
    warning('importPostAcqAnnotations:SchemaMismatch', ...
        ['%s has annotation schema version %d but this GIAnT expects %d; ' ...
         'post-acquisition annotations not used. Align the GIAnT versions used ' ...
         'at the rig and for analysis.'], fnPostAcq, schemaVersion, annotationsSchemaVersion());
    return
end

for DMDix = 1:nDMDs
    if DMDix > numel(postAcq) || isempty(postAcq(DMDix).IM)
        continue %this path was not annotated post-acquisition
    end
    [IMtarget, fn] = getAnnotationImage(trialTable, keepTrials, mocodr, DMDix);
    ROIs(DMDix).fn = fn;
    if isempty(IMtarget)
        warning('importPostAcqAnnotations:NoValidTrials', ...
            'Imaging path %d has no valid trials; post-acquisition annotations not used.', DMDix);
        continue
    end
    [aligned, ok, shiftRC, ncc] = registerAnnotationsToImage(postAcq(DMDix), IMtarget);
    if ~ok
        warning('importPostAcqAnnotations:RegistrationFailed', ...
            'Could not register post-acquisition annotations for imaging path %d (NCC %.2f).', DMDix, ncc);
        continue
    end
    fprintf('Path%d: post-acquisition annotations shifted by [%d %d] (row, col), NCC %.2f\n', ...
        DMDix, shiftRC(1), shiftRC(2), ncc);
    ROIs(DMDix).roiData = aligned.roiData;
    ROIs(DMDix).alignmentShiftRC = shiftRC;
    ROIs(DMDix).alignmentNCC = ncc;
    imported(DMDix) = true;
end
end


function fn = findPostAcqAnnotationFile(datadr)
fn = '';
if isempty(datadr) || ~exist(datadr, 'dir')
    return
end
list = dir(fullfile(datadr, '**', 'postacq_annotations.h5'));
if isempty(list)
    return
end
if numel(list) > 1
    [~, newest] = max([list.datenum]);
    warning('importPostAcqAnnotations:MultipleFiles', ...
        'Found %d postacq_annotations.h5 files under %s; using the newest.', numel(list), datadr);
    list = list(newest);
end
fn = fullfile(list(1).folder, list(1).name);
end
