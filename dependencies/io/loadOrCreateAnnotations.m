function ROIs = loadOrCreateAnnotations(trialTable, keepTrials, mocodr, allowGUI)
%LOADORCREATEANNOTATIONS Load annotations.h5, creating it if needed.
%
% ROIs = loadOrCreateAnnotations(trialTable, keepTrials, mocodr, allowGUI)
%
% If <savedr>/annotations/annotations.h5 exists it is loaded as-is.
% Otherwise, post-acquisition annotations (postacq_annotations.h5 under
% trialTable.datadr) are registered to the motion-corrected data, and the
% user draws ROIs for any imaging path not covered by them. The result is
% saved to annotations.h5.
%
% allowGUI (default: isInteractiveSession()) false never opens the ROI GUI:
% paths not covered by post-acquisition annotations get no user ROIs, with a
% warning, and annotations.h5 is not written, so a later interactive run can
% still draw them. The GUI is never opened in a non-interactive session
% (-batch, -nodisplay), where it would wait forever.

if nargin < 4 || isempty(allowGUI)
    allowGUI = true;
end
allowGUI = allowGUI && isInteractiveSession();

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
skippedDMDs = false(1, numel(imported));
hROIs = cell(1, numel(imported));
for DMDix = find(~imported)
    [IM, fn] = getAnnotationImage(trialTable, keepTrials, mocodr, DMDix);
    if isempty(IM)
        warning('loadOrCreateAnnotations:NoValidTrials', ...
            'Skipping ROI drawing for imaging path %d: no trials passed file verification.', DMDix);
        continue
    end
    ROIs(DMDix).fn = fn;
    if ~allowGUI
        warning('loadOrCreateAnnotations:NoAnnotations', ...
            ['Imaging path %d has no usable annotations and the ROI GUI is not available; ' ...
             'processing it without user ROIs (no soma ROIs, nothing excluded).'], DMDix);
        skippedDMDs(DMDix) = true;
        continue
    end
    hROIs{DMDix} = drawROIs(sqrt(max(0,IM)), annotationsDr, fn);
    drawnDMDs(DMDix) = true;
end
for DMDix = find(drawnDMDs)
    waitfor(hROIs{DMDix}.hF);
    ROIs(DMDix).roiData = hROIs{DMDix}.roiData;
end

if any(skippedDMDs)
    warning('loadOrCreateAnnotations:NotSaved', ...
        ['annotations.h5 not written because imaging path(s) %s were not annotated. ' ...
         'Run annotateROIs or SILo interactively to annotate them.'], mat2str(find(skippedDMDs)));
    return
end
saveAnnotationsH5(fnAnnH5, ROIs);
end
