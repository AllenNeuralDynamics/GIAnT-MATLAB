function [IM, pathIdx] = loadAcquisitionAnnotationImage(fn, maxFrames)
%LOADACQUISITIONANNOTATIONIMAGE Mean image of a raw SLAP2 acquisition tiff for annotation.
%
% [IM, pathIdx] = loadAcquisitionAnnotationImage(fn, maxFrames)
%
% fn is a tiff written by SLAP2 during acquisition (e.g. *_DMD1-TRIAL000001.tif).
% IM is H x W x numChannels in the same orientation as the motion-corrected
% outputs, with unimaged pixels set to nan. pathIdx is the acquisition path
% (DMD) index from the tiff metadata, or nan if absent. At most maxFrames
% (default 2000) evenly spaced frames per channel are averaged. Only those
% pages are read, a chunk at a time, so time and memory do not grow with the
% recording length.

if nargin < 2 || isempty(maxFrames)
    maxFrames = 2000;
end

tiffWarn = warning('off', 'imageio:tiffutils:libtiffWarning'); %some SLAP2 tiffs have an unterminated Software tag
t = Tiff(fn, 'r');
desc = t.getTag('ImageDescription');
t.close();
warning(tiffWarn);
meta = jsondecode(desc);
nchan = meta.numChannels;
if isfield(meta, 'AcquisitionPathIdx')
    pathIdx = double(meta.AcquisitionPathIdx);
else
    pathIdx = nan;
end

try
    [total, imaged, nKept] = accumulateWithIndexedReader(fn, nchan, maxFrames);
catch err
    warning('loadAcquisitionAnnotationImage:IndexedReadFailed', ...
        'Indexed read of %s failed (%s); reading page by page with Tiff instead.', fn, err.message);
    [total, imaged, nKept] = accumulateWithTiff(fn, nchan, maxFrames);
end
IM = single(total / nKept); %plain mean: SLAP2 samples are mostly 0 photons, and trimming would discard the signal
IM(repmat(~imaged, 1, 1, nchan)) = nan;
end


function [total, imaged, nKept] = accumulateWithIndexedReader(fn, nchan, maxFrames)
%ScanImageTiffReader indexes the pages once, then reads just the selected ones
chunkPages = 200;
r = ScanImageTiffReader(fn);
closeReader = onCleanup(@() r.close());
dims = r.shape().dims;
nPages = dims(end);
[pages, nKept] = selectedPages(nPages, nchan, maxFrames);
total = [];
for c0 = 1:chunkPages:numel(pages)
    chunk = pages(c0:min(c0 + chunkPages - 1, numel(pages)));
    data = r.dataPlanes(chunk);
    if isempty(total)
        total = zeros(size(data,1), size(data,2), nchan);
        imaged = false(size(data,1), size(data,2));
    end
    for ch = 1:nchan
        onCh = mod(chunk - 1, nchan) + 1 == ch;
        total(:,:,ch) = total(:,:,ch) + sum(data(:,:,onCh), 3, 'double');
    end
    imaged = imaged | any(data ~= 0, 3); %SLAP2 writes zeros for pixels outside the imaged ROIs
end
total = permute(total, [2 1 3]); %reader returns W x H; registered outputs are H x W
imaged = imaged';
end


function [total, imaged, nKept] = accumulateWithTiff(fn, nchan, maxFrames)
tiffWarn = warning('off', 'imageio:tiffutils:libtiffWarning');
restoreWarn = onCleanup(@() warning(tiffWarn));
t = Tiff(fn, 'r');
closeTiff = onCleanup(@() t.close());
nPages = 1;
while ~t.lastDirectory()
    t.nextDirectory();
    nPages = nPages + 1;
end
[pages, nKept] = selectedPages(nPages, nchan, maxFrames);
wanted = false(1, nPages);
wanted(pages) = true;
t.setDirectory(1);
sz = [t.getTag('ImageLength') t.getTag('ImageWidth')];
total = zeros([sz nchan]);
imaged = false(sz);
for p = 1:pages(end)
    if wanted(p)
        page = t.read();
        ch = mod(p - 1, nchan) + 1;
        total(:,:,ch) = total(:,:,ch) + double(page);
        imaged = imaged | page ~= 0;
    end
    if p < nPages
        t.nextDirectory();
    end
end
end


function [pages, nKept] = selectedPages(nPages, nchan, maxFrames)
%pages of evenly spaced frames; channels are interleaved, one page per channel
nframes = floor(nPages / nchan);
keep = unique(round(linspace(1, nframes, min(nframes, maxFrames))));
pages = reshape((keep(:)' - 1) * nchan + (1:nchan)', 1, []);
nKept = numel(keep);
end
