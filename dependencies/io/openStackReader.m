function [r, planes, nParallel, info, cleanup] = openStackReader(fn, selection, opts)
%OPENSTACKREADER Open a ScanImageTiffReader on fn, chosen for how fn's storage reads best.
%
% [r, planes, nParallel, info, cleanup] = openStackReader(fn, selection, opts)
%
% selection - [] for all planes, a vector of plane indices, or a function
%             handle @(nPlanes) returning plane indices (for callers that need
%             the plane count to choose them)
% opts      - struct with optional fields:
%             Strategy    'auto' (default), 'direct', 'copy', 'pcopy' or 'parallel'
%             ScratchDir  folder for temporary copies (default: GIANT_SCRATCH
%                         environment variable, else tempdir)
%             NumParallel concurrent reads for 'parallel' and 'pcopy' (default 8)
%
% Read the planes with r.dataPlanes(planes, nParallel), or r.data() when all
% planes are wanted and nParallel is 1. Keep cleanup until reading is done:
% clearing it closes the reader and deletes any temporary copy.
%
% Strategy 'auto' follows measurements on SLAP2 recordings:
%   local disk  direct reads
%   SMB share   copy to local scratch first ('copy', the OS copy is ~4x faster
%               than reading over SMB) when >= 25% of the file is needed,
%               otherwise parallel direct reads
%   FUSE (S3)   parallel copy to local scratch ('pcopy'): S3 mounts give ~45
%               MB/s per read stream but ~300 MB/s with 8 concurrent reads, and
%               opening the tiff in place walks ~4,000 page headers one request
%               at a time (~11 s); after the copy that walk is local (<1 s)
%   NFS         direct reads (Code Ocean EFS: ~800 MB/s single stream)
% Without scratch space, FUSE files are read in place with parallel reads.
% On SMB and FUSE storage the page index is predicted from ScanImage's fixed page
% spacing where possible (ScanImageTiffReader 'predict', verified).
% info reports the storage class, strategy and index method used.

if nargin < 2
    selection = [];
end
if nargin < 3 || isempty(opts)
    opts = struct();
end
strategy = lower(getOpt(opts, 'Strategy', 'auto'));
scratchDir = getOpt(opts, 'ScratchDir', '');
if isempty(scratchDir)
    scratchDir = getenv('GIANT_SCRATCH');
end
if isempty(scratchDir)
    scratchDir = tempdir;
end
numParallel = getOpt(opts, 'NumParallel', 8);

info = struct('file', fn, 'storage', storageClass(fn), 'strategy', '', 'copy', '', 'index', '');
if strcmp(strategy, 'auto') && strcmp(info.storage, 'fuse')
    strategy = 'pcopy';
end

r = [];
if strcmp(strategy, 'pcopy') %decided before opening fn: opening it in place is the slow part
    [r, cleanup, info.copy] = openLocalCopy(fn, scratchDir, @(src, dst) parallelCopy(src, dst, numParallel));
    if isempty(r)
        strategy = 'parallel';
    end
end

if isempty(r)
    if any(strcmp(info.storage, {'smb', 'fuse'})) %every IFD read waits on the network
        r = ScanImageTiffReader(fn, 'predict');
    else %local disk, NFS: walking the IFDs is already fast
        r = ScanImageTiffReader(fn, 'walk');
    end
    cleanup = onCleanup(@() r.close());
end
planes = selectPlanes(selection, r.shape().dims(end));

if strcmp(strategy, 'auto')
    switch info.storage
        case 'smb'
            if numel(planes) >= 0.25 * r.shape().dims(end)
                strategy = 'copy';
            else
                strategy = 'parallel';
            end
        otherwise
            strategy = 'direct';
    end
end

if strcmp(strategy, 'copy')
    [rLocal, cleanupLocal, info.copy] = openLocalCopy(fn, scratchDir, @copyWithRetries);
    if isempty(rLocal)
        strategy = 'parallel';
    else
        clear cleanup %closes the reader on the original
        r = rLocal;
        cleanup = cleanupLocal;
    end
end

if strcmp(strategy, 'parallel')
    nParallel = numParallel;
else
    nParallel = 1;
end
info.strategy = strategy;
info.index = r.indexMethod();
end


function [r, cleanup, localCopy] = openLocalCopy(fn, scratchDir, copyFcn)
%copy fn to scratch and open the copy; r is empty (with a warning) if that fails
r = [];
cleanup = [];
localCopy = '';
d = dir(fn);
if freeBytes(scratchDir) < 1.1 * d.bytes + 1e9
    warning('openStackReader:NoScratchSpace', ...
        'Not enough free space in %s to copy %s; reading it in place.', scratchDir, fn);
    return
end
[~, ~, ext] = fileparts(fn);
tmp = [tempname(scratchDir) ext];
try
    copyFcn(fn, tmp);
    r = ScanImageTiffReader(tmp, 'walk');
    cleanup = onCleanup(@() closeAndDelete(r, tmp));
    localCopy = tmp;
catch err
    deleteIfExists(tmp);
    warning('openStackReader:CopyFailed', ...
        'Could not copy %s to %s (%s); reading it in place.', fn, scratchDir, err.message);
end
end


function parallelCopy(src, dst, nParallel)
%copy src to dst with nParallel concurrent readers, each streaming one
%contiguous part of the file; much faster than a single stream on S3 mounts
d = dir(src);
nBytes = d.bytes;
pool = backgroundPool;
nParallel = min(nParallel, pool.NumWorkers);
if nParallel < 2 || nBytes < 64e6
    copyWithRetries(src, dst);
    return
end
raf = java.io.RandomAccessFile(dst, 'rw'); %preallocate (MATLAB's fseek cannot move past
raf.setLength(nBytes);                     %the end of a file), so each part is written in place
raf.close();
edges = round(linspace(0, nBytes, nParallel + 1));
futures(nParallel) = parallel.FevalFuture;
for k = 1:nParallel
    futures(k) = parfeval(pool, @copyRange, 0, src, dst, edges(k), edges(k + 1) - edges(k));
end
wait(futures);
errs = {futures.Error};
failed = find(~cellfun(@isempty, errs), 1);
if ~isempty(failed)
    error('openStackReader:ParallelCopyFailed', 'Part %d failed: %s', failed, errs{failed}.message);
end
d2 = dir(dst);
assert(d2.bytes == nBytes, 'openStackReader:CopyIncomplete', 'Copy is %d of %d bytes', d2.bytes, nBytes);
end


function copyRange(src, dst, offset, nBytes)
fi = fopen(src, 'r');
fo = fopen(dst, 'r+');
closeBoth = onCleanup(@() cellfun(@(f) fclose(f), {fi, fo}));
assert(fi > 0 && fo > 0, 'openStackReader:OpenFailed', 'Could not open %s or %s', src, dst);
fseek(fi, offset, 'bof');
fseek(fo, offset, 'bof');
remaining = nBytes;
while remaining > 0
    buf = fread(fi, min(remaining, 64e6), '*uint8');
    if isempty(buf)
        error('openStackReader:ShortRead', 'Unexpected end of %s at offset %d', src, offset + nBytes - remaining);
    end
    fwrite(fo, buf, 'uint8');
    remaining = remaining - numel(buf);
end
end


function planes = selectPlanes(selection, nPlanes)
if isa(selection, 'function_handle')
    planes = selection(nPlanes);
elseif isempty(selection)
    planes = 1:nPlanes;
else
    planes = selection;
end
end


function v = getOpt(opts, name, default)
if isfield(opts, name) && ~isempty(opts.(name))
    v = opts.(name);
else
    v = default;
end
end


function copyWithRetries(src, dst)
%network copies occasionally fail with transient "system resource" errors
for attempt = 1:5
    try
        copyfile(src, dst);
        return
    catch err
        if attempt == 5
            rethrow(err);
        end
        pause(2 * attempt);
    end
end
end


function b = freeBytes(folder)
b = inf;
try
    if ispc
        d = System.IO.DriveInfo(folder(1:3));
        b = double(d.AvailableFreeSpace);
    else
        [status, out] = system(sprintf('df -B1 --output=avail "%s" | tail -1', folder));
        if status == 0
            b = str2double(strtrim(out));
        end
    end
catch
end
end


function closeAndDelete(r, f)
r.close();
deleteIfExists(f);
end


function deleteIfExists(f)
if exist(f, 'file')
    delete(f);
end
end
