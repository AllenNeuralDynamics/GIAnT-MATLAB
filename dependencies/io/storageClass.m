function cls = storageClass(fn)
%STORAGECLASS Kind of storage a file lives on, for choosing how to read it.
%
% cls = storageClass(fn) returns one of:
%   'local'   local disk
%   'smb'     Windows network share (UNC path or a Network drive), or a CIFS/SMB mount
%   'fuse'    FUSE mount, e.g. S3-backed Code Ocean external data assets
%   'nfs'     NFS mount, e.g. Code Ocean internal data assets (EFS)
%
% Unknown or undetectable storage is treated as 'local'. Set the environment
% variable GIANT_STORAGE_CLASS to one of these values to override detection.
% Results are cached per folder.

override = getenv('GIANT_STORAGE_CLASS');
if ~isempty(override)
    cls = lower(override);
    return
end

persistent cache
if isempty(cache)
    cache = containers.Map('KeyType', 'char', 'ValueType', 'char');
end
folder = fileparts(char(fn));
if isempty(folder)
    folder = pwd;
end
if isKey(cache, folder)
    cls = cache(folder);
    return
end

if ispc
    cls = windowsClass(folder);
else
    cls = unixClass(folder);
end
cache(folder) = cls;
end


function cls = windowsClass(folder)
cls = 'local';
if startsWith(folder, '\\') || startsWith(folder, '//')
    cls = 'smb';
    return
end
try
    drive = System.IO.DriveInfo(folder(1:min(3, end)));
    if strcmp(char(drive.DriveType.ToString()), 'Network')
        cls = 'smb';
    end
catch
    %no .NET, or not a drive path: assume local
end
end


function cls = unixClass(folder)
cls = 'local';
fstype = mountType(folder);
if isempty(fstype)
    return
end
if startsWith(fstype, 'fuse')
    cls = 'fuse';
elseif startsWith(fstype, 'nfs')
    cls = 'nfs';
elseif any(strcmp(fstype, {'cifs', 'smb3', 'smbfs'}))
    cls = 'smb';
end
end


function fstype = mountType(folder)
%file system type of the mount containing folder: the longest mount point
%in /proc/mounts that prefixes it (Linux; no shell needed), else df -T
fstype = '';
folder = char(folder);
try
    target = folder;
    try
        target = char(java.io.File(folder).getCanonicalPath()); %resolve symlinks
    catch
    end
    lines = splitlines(fileread('/proc/mounts'));
    best = 0;
    for k = 1:numel(lines)
        parts = strsplit(strtrim(lines{k}));
        if numel(parts) < 3
            continue
        end
        mnt = strrep(parts{2}, '\040', ' ');
        if (strcmp(target, mnt) || startsWith(target, [regexprep(mnt, '/$', '') '/']) || strcmp(mnt, '/')) ...
                && numel(mnt) > best
            best = numel(mnt);
            fstype = lower(parts{3});
        end
    end
catch
    [status, out] = system(sprintf('df -T "%s" 2>/dev/null | tail -1 | awk ''{print $2}''', folder));
    if status == 0
        fstype = lower(strtrim(out));
    end
end
end
