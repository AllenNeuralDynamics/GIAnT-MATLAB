function fnAnn = postAcquisitionAnnotation(dr, fns)
%POSTACQUISITIONANNOTATION Annotate ROIs (e.g. somata) on raw SLAP2 acquisition files.
%
% fnAnn = postAcquisitionAnnotation(dr)
% fnAnn = postAcquisitionAnnotation(dr, fns)
%
% Intended to run at the rig right after acquisition, while the experimenter
% knows which cells were targeted. Opens the drawROIs GUI on the mean image of
% one acquisition tiff per DMD and saves the ROIs, together with the image they
% were drawn on, to postacq_annotations.h5 next to the selected tiffs.
% annotateROIs/SILo later register these ROIs to the motion-corrected data
% and save them to annotations.h5 (see importPostAcqAnnotations).
%
% dr  - folder to search for acquisition tiffs (default: pwd)
% fns - optional cell array of tiff paths, one per DMD, to skip file selection
%
% Returns the path of the saved file, or '' if the user cancelled.

if nargin < 1 || isempty(dr)
    dr = pwd;
end
dr = char(dr);
fnAnn = '';

if nargin < 2 || isempty(fns)
    fns = selectAcquisitionFiles(dr);
    if isempty(fns)
        disp('Post-acquisition annotation cancelled.')
        return
    end
end
fns = cellstr(fns);

saveDr = fileparts(fns{1});
fnAnnCandidate = fullfile(saveDr, 'postacq_annotations.h5');
if exist(fnAnnCandidate, 'file')
    if ~strcmpi(questdlg('Post-acquisition annotations already exist. Overwrite?', ...
            'Post-acquisition annotation', 'Yes', 'No', 'No'), 'Yes')
        return
    end
end

ROIs = repmat(struct('dr', '', 'fn', '', 'IM', [], 'roiData', {{}}), 1, 0);
hROIs = {};
for fileIx = 1:numel(fns)
    [fileDr, fn, ext] = fileparts(fns{fileIx});
    disp(['Loading ' fn ext ' for annotation...'])
    [IM, pathIdx] = loadAcquisitionAnnotationImage(fns{fileIx});
    if isnan(pathIdx)
        pathIdx = fileIx;
    end
    switch size(IM,3)
        case 2
            IMdisp = IM(:,:,[2 1]); %channels are GR; drawROIs shows ch1 red, ch2 cyan
        otherwise
            IMdisp = IM(:,:,1);
    end
    hROIs{pathIdx} = drawROIs(sqrt(max(0, IMdisp)), fileDr, fn);
    ROIs(pathIdx).dr = fileDr;
    ROIs(pathIdx).fn = [fn ext];
    ROIs(pathIdx).IM = IM;
end
for pathIdx = 1:numel(ROIs)
    if ~isempty(ROIs(pathIdx).IM)
        waitfor(hROIs{pathIdx}.hF);
        ROIs(pathIdx).roiData = hROIs{pathIdx}.roiData;
    end
end

saveAnnotationsH5(fnAnnCandidate, ROIs);
fnAnn = fnAnnCandidate;
disp('Post-acquisition ROI annotations saved to:');
disp(fnAnn);
end


function fns = selectAcquisitionFiles(dr)
%one tiff per DMD; DMD2 is matched to the DMD1 selection by name when possible
fns = {};
default = newestAcquisitionTiff(dr, 1);
[fn1, dr1] = uigetfile(fullfile(fileparts(default), '*DMD1*.tif'), ...
    'Select a DMD1 acquisition file for ANNOTATION', default);
if isequal(fn1, 0)
    return
end
fns{1} = fullfile(dr1, fn1);

fn2 = fullfile(dr1, regexprep(fn1, 'DMD1', 'DMD2', 'once'));
if ~exist(fn2, 'file')
    [fn2, dr2] = uigetfile(fullfile(dr1, '*DMD2*.tif'), ...
        'Select a DMD2 acquisition file for ANNOTATION (cancel to annotate DMD1 only)');
    if isequal(fn2, 0)
        return
    end
    fn2 = fullfile(dr2, fn2);
end
fns{2} = fn2;
end


function fn = newestAcquisitionTiff(dr, DMDix)
%most recent acquisition tiff for this DMD, skipping reference stacks and pipeline outputs
fn = [dr filesep];
list = dir(fullfile(dr, '**', sprintf('*DMD%d*.tif', DMDix)));
list = list(~contains({list.name}, {'REFERENCE', 'REGISTERED'}, 'IgnoreCase', true) ...
    & ~startsWith({list.name}, 'refStack', 'IgnoreCase', true));
if ~isempty(list)
    [~, newest] = max([list.datenum]);
    fn = fullfile(list(newest).folder, list(newest).name);
end
end
