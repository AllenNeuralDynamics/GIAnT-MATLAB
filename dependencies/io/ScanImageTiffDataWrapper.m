function stack = ScanImageTiffDataWrapper(A, tifFile)
%SCANIMAGETIFFDATAWRAPPER Read TIFF stack via ScanImageTiffReader, with fallback.
%
% stack = ScanImageTiffDataWrapper(A, tifFile)
%
% Inputs
%   A       - ScanImageTiffReader object
%   tifFile - path to TIFF file used to construct A
%
% Behavior
%   1) Tries:
%        stack = A.data();
%   2) If that fails, reads the TIFF stack using MATLAB built-ins
%      (readTiffStackImread).
%
% Prefer readStack, which also chooses how to read for the file's storage.

    try
        stack = A.data();
        return
    catch ME
        warning('ScanImageTiffDataWrapper:Fallback', ...
            'A.data() failed (%s). Falling back to built-in TIFF reader.', ...
            ME.message);
    end

    stack = readTiffStackImread(tifFile);
end
