classdef Slap2DataFile < dynamicprops
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
    %     lineLevel      - true: also write the recording as SLAP2 line data
    %                      (hMultiDataFiles, see makeLineData), with known
    %                      tissue motion and a simulated online motion
    %                      correction, for code that reads lines and headers
    %     linesPerCycle, lineMotion, updateEvery - see makeLineData

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
        truthMotion = []      % lineLevel: tissue motion [rows cols] at every line
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
            if obj.config.lineLevel
                % a dynamic property, so that only lineLevel fakes look like a
                % reader with line data to isprop(..., 'hMultiDataFiles')
                addprop(obj, 'hMultiDataFiles');
                [obj.hMultiDataFiles, obj.truthMotion] = slap2.Slap2DataFile.makeLineData( ...
                    filename, obj.config, obj.texture, obj.totalNumLines);
            end
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
                'firstLineTime', NaT, 'lineLevel', false, 'linesPerCycle', 20, ...
                'lineMotion', struct('drift', [3 -2], 'amp', [1.8 1.5], 'period', [3200 2700]), ...
                'updateEvery', 26);
            sidecar = [filename '.fake.mat'];
            if isfile(sidecar)
                s = load(sidecar, 'cfg');
                for f = fieldnames(s.cfg)'
                    cfg.(f{1}) = s.cfg.(f{1});
                end
            end
        end

        function [M, m] = makeLineData(filename, cfg, tex, nLines)
            % A SLAP2 data file for the recording, laid out as SLAP2 writes it
            % (file header, then cycles of lines: a 40-byte line header with
            % flags, xOffset_pix, yOffset_pix and zOffset_um, then the samples
            % of each channel), and the hMultiDataFiles fields that the readers
            % use. Each of the linesPerCycle lines samples every
            % linesPerCycle-th ROI pixel. The tissue moves by
            %   m(l) = drift*l/nLines + amp.*sin(2*pi*l./period)   [rows cols]
            % and a simulated online correction sets the header offsets every
            % updateEvery lines to x = round(m_cols), y = round(m_rows). As
            % measured on SLAP2, a new y reaches the DMD 2 lines before the
            % header records it (6 when x changes too). A sample of DMD pixel
            % (r, c) on line l is the texture at (r + yDmd - m_rows, c + x -
            % m_cols), times 10, divided by the channel number.
            P = cfg.linesPerCycle; nC = nLines / P;
            assert(nC == round(nC), 'lineLevel fakes need a whole number of cycles.');
            NC = cfg.frameSize(1); NR = cfg.frameSize(2); nCh = cfg.numChannels;
            [cc, rr] = meshgrid(cfg.roiCols, cfg.roiRows);
            ids = (rr(:) - 1) * NR + cc(:);                 % 1-based raster IDs of the ROI
            spIDs = arrayfun(@(L) uint32(ids(L:P:end) - 1), 1:P, 'UniformOutput', false);
            l = (1:nLines)';
            mo = cfg.lineMotion;
            m = l / nLines .* mo.drift + mo.amp .* sin(2 * pi * l ./ mo.period);
            x = zeros(nLines, 1); y = zeros(nLines, 1);
            for s = 1:cfg.updateEvery:nLines
                x(s:end) = round(m(s, 2)); y(s:end) = round(m(s, 1));
            end
            yDmd = y;
            for s = (find(diff(y) ~= 0) + 1)'
                k = 2; if x(s) ~= x(s-1), k = 6; end
                yDmd(max(1, s - k):s-1) = y(s);
            end
            % file layout
            H = 10; hw = 20; nSp = cellfun(@numel, spIDs(:));
            lineWords = hw + nCh * nSp;
            lineHeaderIdxs = H + 1 + [0; cumsum(lineWords(1:end-1))];
            cycleWords = sum(lineWords);
            words = zeros(H + cycleWords * nC, 1, 'int16');
            [tc, tr] = meshgrid(1:size(tex, 2), 1:size(tex, 1));
            for c = 1:nC
                for L = 1:P
                    li = (c - 1) * P + L;
                    at = (c - 1) * cycleWords + lineHeaderIdxs(L);
                    words(at + (10:11)) = int16([x(li) y(li)]);
                    id = double(spIDs{L}) + 1;
                    r = floor((id - 1) / NR) + 1; col = mod(id - 1, NR) + 1;
                    v = interp2(tc, tr, tex, col + x(li) - m(li, 2), r + yDmd(li) - m(li, 1), 'linear', 10);
                    for ch = 1:nCh
                        words(at + hw + (ch - 1) * nSp(L) + (0:nSp(L)-1)) = int16(round(10 * v / ch));
                    end
                end
            end
            datFile = [filename '.lines.dat'];
            fid = fopen(datFile, 'w'); fwrite(fid, words, 'int16'); fclose(fid);
            header = struct('linesPerCycle', P, 'dmdPixelsPerRow', NR, 'dmdPixelsPerColumn', NC, ...
                'numChannels', nCh, 'channelsInterleave', false, 'bytesPerCycle', 2 * cycleWords, ...
                'lineHeaderSizeBytes', 2 * hw);
            D = struct('datFileName', datFile, 'rawData', struct('Format', 'int16', 'Data', words), ...
                'header', header, 'lineHeaderIdxs', lineHeaderIdxs, 'numCycles', nC);
            M = struct('header', header, 'filename', datFile, 'totalNumLines', nLines, 'numCycles', nC, ...
                'fastZs', 0, 'lineFastZIdxs', ones(1, P), 'lineSuperPixelIDs', {spIDs}, ...
                'lineHeaderIdxs', lineHeaderIdxs, 'lineDataStartIdxs', lineHeaderIdxs + hw, ...
                'lineDataNumElements', nCh * nSp, 'hDataFiles', D);
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
