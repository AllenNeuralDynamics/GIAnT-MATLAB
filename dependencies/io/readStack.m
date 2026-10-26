function [data, meta, desc, info] = readStack(fn, planes, varargin)
%READSTACK Read an imaging stack (ScanImage tiff or h5 movie), efficiently for its storage.
%
% data = readStack(fn)          whole stack
% data = readStack(fn, planes)  only the given planes (tiff pages, or indices
%                               along the last dimension of an h5 dataset)
% [data, meta, desc, info] = readStack(fn, planes, 'Name', value, ...)
%
% Tiffs: data has the layout and class of ScanImageTiffReader.data() (native
% type, W x H x planes); meta is the reader's metadata and desc the per-page
% ImageDescription strings (read only when requested). Compressed tiffs, or
% files the reader cannot parse, are read with imread instead.
% h5 files: data is as from h5read (first dataset unless 'Dataset' is given);
% meta and desc are empty.
%
% Options: 'Strategy' ('auto', 'direct', 'copy', 'parallel'), 'ScratchDir',
% 'NumParallel', 'Dataset' (h5). How 'auto' chooses per storage (local, SMB,
% S3/FUSE, NFS) is described in openStackReader. info reports the storage
% class and strategy used.

if nargin < 2
    planes = [];
end
opts = parseOptions(varargin);
meta = '';
desc = {};

[~, ~, ext] = fileparts(fn);
if any(strcmpi(ext, {'.h5', '.hdf5'}))
    [data, info] = readH5(fn, planes, opts);
    return
end

try
    [r, sel, nParallel, info, cleanup] = openStackReader(fn, planes, opts); %#ok<ASGLU> cleanup closes r
    if r.compression() ~= 1
        error('readStack:Compressed', 'Compressed tiff (compression %d).', r.compression());
    end
    if isempty(planes) && nParallel == 1
        data = r.data();
    else
        data = r.dataPlanes(sel, nParallel);
    end
    meta = r.metadata();
    if nargout > 2
        desc = r.descriptions();
    end
catch err
    warning('readStack:Fallback', 'Reading %s with imread (%s).', fn, err.message);
    [data, desc] = readTiffStackImread(fn, planes);
    info = struct('file', fn, 'storage', storageClass(fn), 'strategy', 'imread', 'copy', '', 'index', '');
end
end


function [data, info] = readH5(fn, planes, opts)
info = struct('file', fn, 'storage', storageClass(fn), 'strategy', 'direct', 'copy', '', 'index', '');
ds = opts.Dataset;
if isempty(ds)
    h = h5info(fn);
    ds = h.Datasets(1).Name;
end
ds = ['/' regexprep(ds, '^/', '')];
if isempty(planes)
    data = h5read(fn, ds);
    return
end
sz = h5info(fn, ds).Dataspace.Size;
lo = min(planes);
start = ones(1, numel(sz)); start(end) = lo;
count = sz; count(end) = max(planes) - lo + 1;
block = h5read(fn, ds, start, count); %one hyperslab covering the requested planes
idx = repmat({':'}, 1, numel(sz));
idx{end} = planes - lo + 1;
data = block(idx{:});
end


function opts = parseOptions(args)
p = inputParser;
p.addParameter('Strategy', 'auto');
p.addParameter('ScratchDir', '');
p.addParameter('NumParallel', 8);
p.addParameter('Dataset', '');
p.parse(args{:});
opts = p.Results;
end
