function ROIs = loadOrCreateAnnotations(trialTable, keepTrials, mocodr, allowGUI)
%LOADORCREATEANNOTATIONS Load annotations.h5, creating or completing it if needed.
%
% ROIs = loadOrCreateAnnotations(trialTable, keepTrials, mocodr, allowGUI)
%
% Imaging paths are annotated from, in order: <savedr>/annotations/annotations.h5;
% post-acquisition annotations (postacq_annotations.h5 under trialTable.datadr),
% registered to the motion-corrected data; or ROIs the user draws. The result
% is saved to annotations.h5. A path entry with an empty fn is "not annotated
% yet" (e.g. it had no valid trials, or the GUI was unavailable); such paths
% are annotated again on later calls once possible, and the rest of the file
% is kept as is.
%
% allowGUI (default true) false never opens the ROI GUI; the GUI is also never
% opened in a non-interactive session (-batch, -nodisplay), where it would wait
% forever. Paths that would need it get no user ROIs, with a warning.

if nargin < 4 || isempty(allowGUI)
    allowGUI = true;
end
allowGUI = allowGUI && isInteractiveSession();

annotationsDr = fullfile(trialTable.savedr, 'annotations');
if ~exist(annotationsDr, 'dir')
    mkdir(annotationsDr);
end
fnAnnH5 = fullfile(annotationsDr, 'annotations.h5');
nDMDs = size(trialTable.filename, 1);
hasValidTrials = any(keepTrials, 2)';

saved = [];
if exist(fnAnnH5, 'file')
    saved = loadAnnotationsH5(fnAnnH5);
    annotated = false(1, nDMDs);
    for DMDix = 1:min(nDMDs, numel(saved))
        annotated(DMDix) = ~isempty(saved(DMDix).fn);
    end
    if ~any(~annotated & hasValidTrials)
        ROIs = saved;
        return
    end
end

[ROIs, resolved] = importPostAcqAnnotations(trialTable, keepTrials, mocodr);
for DMDix = 1:min(nDMDs, numel(saved)) %keep paths already annotated in annotations.h5
    if ~isempty(saved(DMDix).fn)
        for fld = fieldnames(saved(DMDix))' %all of it, including any stored image
            ROIs(DMDix).(fld{1}) = saved(DMDix).(fld{1});
        end
        resolved(DMDix) = true;
    end
end

drawnDMDs = false(1, nDMDs);
skippedDMDs = false(1, nDMDs);
hROIs = cell(1, nDMDs);
for DMDix = find(~resolved)
    [IM, fn] = getAnnotationImage(trialTable, keepTrials, mocodr, DMDix);
    if isempty(IM)
        warning('loadOrCreateAnnotations:NoValidTrials', ...
            'Imaging path %d has no valid trials; it will be annotated once it has.', DMDix);
        markUnannotated(DMDix);
        continue
    end
    if ~allowGUI
        warning('loadOrCreateAnnotations:NoAnnotations', ...
            ['Imaging path %d has no usable annotations and the ROI GUI is not available; ' ...
             'processing it without user ROIs (no soma ROIs, nothing excluded).'], DMDix);
        markUnannotated(DMDix);
        skippedDMDs(DMDix) = true;
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

if any(skippedDMDs)
    warning('loadOrCreateAnnotations:Unannotated', ...
        ['Imaging path(s) %s saved to annotations.h5 as not annotated; run annotateROIs or ' ...
         'SILo interactively to annotate them.'], mat2str(find(skippedDMDs)));
end
saveAnnotationsH5(fnAnnH5, ROIs);

    function markUnannotated(DMDix)
        ROIs(DMDix).fn = ''; %empty fn marks the path for annotation on a later call
        ROIs(DMDix).roiData = {};
        ROIs(DMDix).alignmentShiftRC = [];
        ROIs(DMDix).alignmentNCC = [];
    end
end
