function ROIs = loadOrCreateAnnotations(trialTable, keepTrials, mocodr)
%LOADORCREATEANNOTATIONS Load annotations.h5, creating it if needed.
%
% ROIs = loadOrCreateAnnotations(trialTable, keepTrials, mocodr)
%
% If <savedr>/annotations/annotations.h5 exists it is loaded as-is.
% Otherwise, post-acquisition annotations (postacq_annotations.h5 under
% trialTable.datadr) are registered to the motion-corrected data, the user
% draws ROIs for any imaging path not covered by them, and the result is
% saved to annotations.h5.

annotationsDr = fullfile(trialTable.savedr, 'annotations');
if ~exist(annotationsDr, 'dir')
    mkdir(annotationsDr);
end
fnAnnH5 = fullfile(annotationsDr, 'annotations.h5');
if exist(fnAnnH5, 'file')
    ROIs = loadAnnotationsH5(fnAnnH5);
    return
end

[ROIs, imported] = importPostAcqAnnotations(trialTable, keepTrials, mocodr);

drawnDMDs = false(1, numel(imported));
hROIs = cell(1, numel(imported));
for DMDix = find(~imported)
    [IM, fn] = getAnnotationImage(trialTable, keepTrials, mocodr, DMDix);
    if isempty(IM)
        warning('loadOrCreateAnnotations:NoValidTrials', ...
            'Skipping ROI drawing for imaging path %d: no trials passed file verification.', DMDix);
        continue
    end
    hROIs{DMDix} = drawROIs(sqrt(max(0,IM)), annotationsDr, fn);
    ROIs(DMDix).fn = fn;
    drawnDMDs(DMDix) = true;
end
for DMDix = find(drawnDMDs)
    waitfor(hROIs{DMDix}.hF);
    ROIs(DMDix).roiData = hROIs{DMDix}.roiData;
end
saveAnnotationsH5(fnAnnH5, ROIs);
end
