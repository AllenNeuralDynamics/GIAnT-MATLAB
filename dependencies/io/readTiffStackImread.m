function [stack, desc] = readTiffStackImread(tifFile, planes)
%READTIFFSTACKIMREAD Read tiff pages with MATLAB's imread, in ScanImageTiffReader.data() layout.
%
% [stack, desc] = readTiffStackImread(tifFile, planes)
%
% Fallback for tiffs ScanImageTiffReader cannot decode (e.g. compressed).
% planes ([] = all) selects pages. stack has the reader's layout: W x H x N
% for single-sample pages, C x W x H x N for multi-sample pages, native type.
% desc holds each selected page's ImageDescription ('' if absent).

info = imfinfo(tifFile);
if isempty(info)
    error('readTiffStackImread:EmptyTiff', 'No TIFF frames found in file: %s', tifFile);
end
if nargin < 2 || isempty(planes)
    planes = 1:numel(info);
end

first = imread(tifFile, planes(1), 'Info', info);
if ismatrix(first)
    stack = zeros(size(first, 2), size(first, 1), numel(planes), class(first));
    for k = 1:numel(planes)
        stack(:, :, k) = imread(tifFile, planes(k), 'Info', info)';
    end
else
    stack = zeros(size(first, 3), size(first, 2), size(first, 1), numel(planes), class(first));
    for k = 1:numel(planes)
        stack(:, :, :, k) = permute(imread(tifFile, planes(k), 'Info', info), [3 2 1]);
    end
end

desc = cell(numel(planes), 1);
for k = 1:numel(planes)
    if isfield(info, 'ImageDescription') && ~isempty(info(planes(k)).ImageDescription)
        desc{k} = info(planes(k)).ImageDescription;
    else
        desc{k} = '';
    end
end
end
