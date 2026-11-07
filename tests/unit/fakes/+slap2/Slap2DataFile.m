classdef Slap2DataFile < handle
    %SLAP2DATAFILE Test double for the SLAP2 reader.
    %   Reads the line count from the .dat file's text contents, so tests
    %   can create .dat files containing a single integer.
    %
    %   getImage returns a smooth random texture inside a fixed ROI and NaN
    %   elsewhere. An optional sidecar, [filename '.fake.mat'], holds a
    %   struct cfg that overrides the defaults in loadConfig:
    %     numChannels    - channel c is the texture scaled by 1/c
    %     frameSize      - [rows cols] of each returned image
    %     roiRows/Cols   - pixels inside the ROI
    %     shiftAfterLine - from this line on, sample the texture at +shift
    %     shift          - [rows cols] integer shift
    %     failOnCall     - getImage errors on this call number and later
    %     firstLineTime  - datetime of the file's first line; if set, the
    %                      file gets FPGA line-clock metadata like the real
    %                      reader's (firstLineTimestamp, fpgaTimeReference)

    properties (Constant)
        LinePeriod_s = 1e-4
        FpgaClock_Hz = 2e8
        FpgaTimeReference = datetime(2024, 2, 1, 9, 0, 0)
    end

    properties
        filename
        totalNumLines
        numChannels
        hDataFile
    end

    properties (Access = private)
        config
        texture
        numGetImageCalls = 0
    end

    methods
        function obj = Slap2DataFile(filename)
            obj.filename = filename;
            obj.totalNumLines = str2double(strtrim(fileread(filename)));
            obj.config = slap2.Slap2DataFile.loadConfig(filename);
            obj.numChannels = obj.config.numChannels;
            obj.texture = slap2.Slap2DataFile.makeTexture(obj.config.frameSize);
            linesPerCycle = 1000;
            metaData = struct('linePeriod_s', slap2.Slap2DataFile.LinePeriod_s);
            header = struct();
            firstLineTimestamp = uint64(0);
            if ~isnat(obj.config.firstLineTime)
                metaData.fpgaTimeReference = slap2.Slap2DataFile.FpgaTimeReference;
                header.fpgaSystemClock_Hz = slap2.Slap2DataFile.FpgaClock_Hz;
                firstLineTimestamp = uint64(round(seconds(obj.config.firstLineTime ...
                    - slap2.Slap2DataFile.FpgaTimeReference) * slap2.Slap2DataFile.FpgaClock_Hz));
            end
            obj.hDataFile = struct( ...
                'metaData', metaData, ...
                'header', header, ...
                'firstLineTimestamp', firstLineTimestamp, ...
                'fastZs', 0, ...
                'lineHeaderIdxs', 1:linesPerCycle, ...
                'numCycles', ceil(obj.totalNumLines/linesPerCycle), ...
                'getLineHeader', @(~, ~) struct('xOffset_pix', 0, ...
                    'yOffset_pix', 0, 'zOffset_um', 0));
        end

        function [IM, sampleTimes, freshness] = getImage(obj, channel, lineIdx, ~, ~, ~)
            obj.numGetImageCalls = obj.numGetImageCalls + 1;
            cfg = obj.config;
            if obj.numGetImageCalls >= cfg.failOnCall
                error('Slap2DataFile:FakeReadError', 'Simulated read failure.');
            end
            shift = [0 0];
            if lineIdx >= cfg.shiftAfterLine
                shift = cfg.shift;
            end
            IM = nan(cfg.frameSize);
            IM(cfg.roiRows, cfg.roiCols) = obj.texture( ...
                cfg.roiRows + shift(1), cfg.roiCols + shift(2)) / channel;
            sampleTimes = [];
            freshness = ones(cfg.frameSize);
        end
    end

    methods (Static, Access = private)
        function cfg = loadConfig(filename)
            cfg = struct('numChannels', 1, 'frameSize', [40 40], ...
                'roiRows', 11:30, 'roiCols', 11:30, ...
                'shiftAfterLine', Inf, 'shift', [0 0], 'failOnCall', Inf, ...
                'firstLineTime', NaT);
            sidecar = [filename '.fake.mat'];
            if isfile(sidecar)
                s = load(sidecar, 'cfg');
                for f = fieldnames(s.cfg)'
                    cfg.(f{1}) = s.cfg.(f{1});
                end
            end
        end

        function tex = makeTexture(frameSize)
            % Smooth random texture in [10, 100], identical for every file.
            stream = RandStream('mt19937ar', 'Seed', 0);
            noise = randn(stream, frameSize + 12);
            [x, y] = meshgrid(-6:6);
            kernel = exp(-(x.^2 + y.^2)/8);
            tex = conv2(noise, kernel/sum(kernel(:)), 'valid');
            tex = 10 + 90*(tex - min(tex(:)))/(max(tex(:)) - min(tex(:)));
        end
    end
end
