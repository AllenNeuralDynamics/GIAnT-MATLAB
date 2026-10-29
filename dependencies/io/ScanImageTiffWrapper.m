function [data, meta, desc] = ScanImageTiffWrapper(fn)
%SCANIMAGETIFFWRAPPER Read a whole ScanImage tiff: data/meta/desc as from ScanImageTiffReader.
%   desc - cell array of ImageDescription strings per frame (read only if requested)
%   Reads through readStack, which picks the fastest way to read for the
%   file's storage (local disk, SMB share, S3/FUSE or NFS mount).
if nargout > 2
    [data, meta, desc] = readStack(fn);
else
    [data, meta] = readStack(fn);
end
end
