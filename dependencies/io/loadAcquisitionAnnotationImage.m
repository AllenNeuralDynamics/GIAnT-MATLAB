function [IM, pathIdx] = loadAcquisitionAnnotationImage(fn, maxFrames)
%LOADACQUISITIONANNOTATIONIMAGE Mean image of a raw SLAP2 acquisition tiff for annotation.
%
% [IM, pathIdx] = loadAcquisitionAnnotationImage(fn, maxFrames)
%
% fn is a tiff written by SLAP2 during acquisition (e.g. *_DMD1-TRIAL000001.tif).
% IM is H x W x numChannels in the same orientation as the motion-corrected
% outputs, with unimaged pixels set to nan. pathIdx is the acquisition path
% (DMD) index from the tiff metadata, or nan if absent. At most maxFrames
% (default 2000) evenly spaced frames per channel are averaged.

if nargin < 2 || isempty(maxFrames)
    maxFrames = 2000;
end

[data, ~, desc] = ScanImageTiffWrapper(fn);
meta = jsondecode(desc{1});
nchan = meta.numChannels;
if isfield(meta, 'AcquisitionPathIdx')
    pathIdx = double(meta.AcquisitionPathIdx);
else
    pathIdx = nan;
end

nframes = floor(size(data,3) ./ nchan);
data = reshape(data(:,:,1:nframes*nchan), size(data,1), size(data,2), nchan, nframes);
unimaged = all(data == 0, [3 4]); %SLAP2 writes zeros for pixels outside the imaged ROIs
keep = unique(round(linspace(1, nframes, min(nframes, maxFrames))));
%plain mean, not a trimmed mean: SLAP2 samples are mostly 0 photons, so the
%signal is in the rare nonzero samples that trimming would discard
IM = squeeze(mean(single(data(:,:,:,keep)), 4));
IM(repmat(unimaged, 1, 1, nchan)) = nan;
IM = permute(IM, [2 1 3]); %tiff reader returns W x H; registered outputs are H x W
end
