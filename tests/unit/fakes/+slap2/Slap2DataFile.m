classdef Slap2DataFile < handle
    %SLAP2DATAFILE Test double for the SLAP2 reader.
    %   Reads the line count from the .dat file's text contents, so tests
    %   can create .dat files containing a single integer.

    properties (Constant)
        LinePeriod_s = 1e-4
    end

    properties
        filename
        totalNumLines
        hDataFile
    end

    methods
        function obj = Slap2DataFile(filename)
            obj.filename = filename;
            obj.totalNumLines = str2double(strtrim(fileread(filename)));
            obj.hDataFile = struct('metaData', ...
                struct('linePeriod_s', slap2.Slap2DataFile.LinePeriod_s));
        end
    end
end
