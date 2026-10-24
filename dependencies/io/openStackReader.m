function [r, planes, nParallel, info, cleanup] = openStackReader(fn, selection, opts)
%OPENSTACKREADER Open a ScanImageTiffReader on fn, chosen for how fn's storage reads best.
%
% [r, planes, nParallel, info, cleanup] = openStackReader(fn, selection, opts)
%
% selection - [] for all planes, a vector of plane indices, or a function
%             handle @(nPlanes) returning plane indices (for callers that need
%             the plane count to choose them)
% opts      - struct with optional fields:
%             Strategy    'auto' (default), 'direct', 'copy' or 'parallel'
%             ScratchDir  folder for temporary copies (default: GIANT_SCRATCH
%                         environment variable, else tempdir)
%             NumParallel concurrent reads for 'parallel' (default 8)
%
% Read the planes with r.dataPlanes(planes, nParallel), or r.data() when all
% planes are wanted and nParallel is 1. Keep cleanup until reading is done:
% clearing it closes the reader and deletes any temporary copy.
%
% Strategy 'auto' follows measurements on SLAP2 recordings:
%   local disk     direct reads
%   SMB share      copy to local scratch first when >= 25% of the file is
%                  needed (the OS copy is ~4x faster than reading over SMB),
%                  otherwise parallel direct reads
%   FUSE (S3)      parallel direct reads: S3 mounts give ~40-50 MB/s per read
%                  stream and do not cache, so copying does not help; 8
%                  concurrent reads reach ~170 MB/s on Code Ocean
%   NFS            direct reads (Code Ocean EFS: ~800 MB/s single stream,
%                  slower with concurrent reads)
% info reports the storage class and strategy used.

if nargin < 2
    selection = [];
end
if nargin < 3 || isempty(opts)
    opts = struct();
end
strategy = getOpt(opts, 'Strategy', 'auto');
scratchDir = getOpt(opts, 'ScratchDir', '');
if isempty(scratchDir)
    scratchDir = getenv('GIANT_SCRATCH');
end
if isempty(scratchDir)
    scratchDir = tempdir;
end
numParallel = getOpt(opts, 'NumParallel', 8);

info = struct('file', fn, 'storage', storageClass(fn), 'strategy', '', 'copy', '');
r = ScanImageTiffReader(fn);
cleanup = onCleanup(@() r.close());
nPlanes = r.shape().dims(end);
if isa(selection, 'function_handle')
    planes = selection(nPlanes);
elseif isempty(selection)
    planes = 1:nPlanes;
else
    planes = selection;
end

if strcmpi(strategy, 'auto')
    switch info.storage
        case 'smb'
            if numel(planes) >= 0.25 * nPlanes
                strategy = 'copy';
            else
                strategy = 'parallel';
            end
        case 'fuse'
            strategy = 'parallel';
        otherwise
            strategy = 'direct';
    end
end

if strcmpi(strategy, 'copy')
    d = dir(fn);
    if freeBytes(scratchDir) < 1.1 * d.bytes + 1e9
        warning('openStackReader:NoScratchSpace', ...
            'Not enough free space in %s to copy %s; reading it directly.', scratchDir, fn);
        strategy = 'parallel';
    else
        [~, ~, ext] = fileparts(fn);
        localCopy = [tempname(scratchDir) ext];
        try
            copyWithRetries(fn, localCopy);
            r.close();
            r = ScanImageTiffReader(localCopy);
            cleanup = onCleanup(@() closeAndDelete(r, localCopy));
            info.copy = localCopy;
        catch err
            deleteIfExists(localCopy);
            warning('openStackReader:CopyFailed', ...
                'Could not copy %s to %s (%s); reading it directly.', fn, scratchDir, err.message);
            strategy = 'parallel';
        end
    end
end

switch lower(strategy)
    case 'parallel'
        nParallel = numParallel;
    otherwise
        nParallel = 1;
end
info.strategy = lower(strategy);
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
