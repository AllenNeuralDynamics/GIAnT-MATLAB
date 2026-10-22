function [IM, pathIdx] = loadAcquisitionAnnotationImage(fn, maxFrames, method)
%LOADACQUISITIONANNOTATIONIMAGE Mean image of a raw SLAP2 acquisition tiff for annotation.
%
% [IM, pathIdx] = loadAcquisitionAnnotationImage(fn, maxFrames, method)
%
% fn is a tiff written by SLAP2 during acquisition (e.g. *_DMD1-TRIAL000001.tif).
% IM is H x W x numChannels in the same orientation as the motion-corrected
% outputs, with unimaged pixels set to nan. pathIdx is the acquisition path
% (DMD) index from the tiff metadata, or nan if absent. At most maxFrames
% (default 2000) evenly spaced frames per channel are averaged.
%
% method 'auto' (default) reads the whole file at once when it fits easily in
% available memory (fastest), otherwise one page at a time so memory use does
% not grow with the recording length. 'whole' and 'pages' force either.

if nargin < 2 || isempty(maxFrames)
    maxFrames = 2000;
end
if nargin < 3 || isempty(method)
    method = 'auto';
end
if strcmpi(method, 'auto')
    method = 'pages';
    if fitsInMemory(fn)
        method = 'whole';
    end
end

if strcmpi(method, 'whole')
    [IM, pathIdx] = readWholeFile(fn, maxFrames);
    return
end
try
    [IM, pathIdx] = readSelectedPages(fn, maxFrames);
catch err
    warning('loadAcquisitionAnnotationImage:PageReadFailed', ...
        'Page-wise read of %s failed (%s); reading the whole file instead.', fn, err.message);
    [IM, pathIdx] = readWholeFile(fn, maxFrames);
end
end


function tf = fitsInMemory(fn)
%the whole-file read holds the stack plus a zero mask and a single-precision
%copy of the selected frames: about 4x the file size; keep a 2x margin
tf = false;
if ~ispc
    return %no portable available-memory query; read page-wise
end
d = dir(fn);
try
    m = memory;
    tf = 8 * d.bytes < m.MemAvailableAllArrays;
catch
end
end


function [IM, pathIdx] = readSelectedPages(fn, maxFrames)
t = Tiff(fn, 'r');
closeTiff = onCleanup(@() t.close());
[nchan, pathIdx] = parseMeta(t.getTag('ImageDescription'));

nPages = 1;
while ~t.lastDirectory()
    t.nextDirectory();
    nPages = nPages + 1;
end
nframes = floor(nPages / nchan);
keep = frameSelection(nframes, maxFrames);
pageWanted = false(1, nPages);
pageWanted(reshape((keep(:)' - 1) * nchan + (1:nchan)', 1, [])) = true;

t.setDirectory(1);
sz = [t.getTag('ImageLength') t.getTag('ImageWidth')];
total = zeros([sz nchan]);
imaged = false(sz);
for p = 1:find(pageWanted, 1, 'last')
    if pageWanted(p)
        page = t.read();
        ch = mod(p - 1, nchan) + 1;
        total(:,:,ch) = total(:,:,ch) + double(page);
        imaged = imaged | page ~= 0; %SLAP2 writes zeros for pixels outside the imaged ROIs
    end
    if p < nPages
        t.nextDirectory();
    end
end
IM = single(total / numel(keep)); %plain mean: SLAP2 samples are mostly 0 photons, and trimming would discard the signal
IM(repmat(~imaged, 1, 1, nchan)) = nan;
end


function [IM, pathIdx] = readWholeFile(fn, maxFrames)
[data, ~, desc] = ScanImageTiffWrapper(fn);
[nchan, pathIdx] = parseMeta(desc{1});
nframes = floor(size(data,3) ./ nchan);
data = reshape(data(:,:,1:nframes*nchan), size(data,1), size(data,2), nchan, nframes);
keep = frameSelection(nframes, maxFrames);
data = data(:,:,:,keep);
unimaged = all(data == 0, [3 4]);
IM = squeeze(mean(single(data), 4));
IM(repmat(unimaged, 1, 1, nchan)) = nan;
IM = permute(IM, [2 1 3]); %tiff reader returns W x H; registered outputs are H x W
end


function [nchan, pathIdx] = parseMeta(desc)
meta = jsondecode(desc);
nchan = meta.numChannels;
if isfield(meta, 'AcquisitionPathIdx')
    pathIdx = double(meta.AcquisitionPathIdx);
else
    pathIdx = nan;
end
end


function keep = frameSelection(nframes, maxFrames)
keep = unique(round(linspace(1, nframes, min(nframes, maxFrames))));
end
