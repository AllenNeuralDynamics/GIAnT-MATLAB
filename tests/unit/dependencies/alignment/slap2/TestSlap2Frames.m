classdef TestSlap2Frames < matlab.unittest.TestCase
    %TESTSLAP2FRAMES Unit tests for slap2Frames and slap2FrameSums.
    %   Uses a struct in place of slap2.Slap2DataFile. Its hMultiDataFiles
    %   holds a small parse plan; getLineData returns a deterministic value
    %   for every (line, cycle, superpixel, channel) and, like the real
    %   reader, returns the channels min:max in file order. makeDiskData also
    %   writes these values to two data files laid out as SLAP2 writes them, so
    %   that the file reader is tested too. Expected images come from a
    %   brute-force loop that follows getImageBatched (or getImages).
    %   Run with: runtests('tests/unit/dependencies/alignment/slap2/TestSlap2Frames.m')

    properties (Constant)
        NR = 12   % dmdPixelsPerRow (image columns)
        NC = 10   % dmdPixelsPerColumn (image rows)
        P = 6     % lines per cycle
        nCycles = 20
    end

    methods (TestClassSetup)
        function addPaths(testCase)
            import matlab.unittest.fixtures.PathFixture
            testsDir = fileparts(mfilename('fullpath'));
            repoRoot = fileparts(fileparts(fileparts(fileparts(fileparts(testsDir)))));
            testCase.applyFixture(PathFixture(fullfile(repoRoot, 'dependencies', 'alignment', 'slap2')));
        end
    end

    methods (Test)
        function testMatchesGetImageWeighting(testCase)
            S2 = testCase.makeData();
            times = [1 7 30 61 120];
            [IM, F] = slap2Frames(S2, [1 2], times, 4);
            [IMref, Fref] = testCase.reference(S2, [1 2], times, 4, 'exp', zeros(testCase.totalLines, 2), [0 1; 1 0]);
            testCase.verifyEqual(isnan(IM), isnan(IMref));
            testCase.verifyEqual(double(IM(~isnan(IM))), IMref(~isnan(IMref)), 'RelTol', 1e-6);
            testCase.verifyEqual(double(F), Fref, 'RelTol', 1e-6);
        end

        function testChannelsInRequestedOrder(testCase)
            S2 = testCase.makeData();
            IM12 = slap2Frames(S2, [1 2], 40, 4);
            IM21 = slap2Frames(S2, [2 1], 40, 4);
            testCase.verifyEqual(IM21(:, :, 1), IM12(:, :, 2));
            testCase.verifyEqual(IM21(:, :, 2), IM12(:, :, 1));
            IM2 = slap2Frames(S2, 2, 40, 4);
            testCase.verifyEqual(IM2, IM12(:, :, 2));
        end

        function testConstantOffsetIsIntegerShift(testCase)
            S2 = testCase.makeData();
            off = testCase.offsets(2 * ones(testCase.totalLines, 1), -1 * ones(testCase.totalLines, 1));
            plain = slap2Frames(S2, 1, 50, 4);
            moved = slap2Frames(S2, 1, 50, 4, CompensateOnline = true, Offsets = off, RefOffset = [0 0]);
            % default mapping: rows follow y (-1), columns follow x (+2)
            expected = nan(size(plain), 'single');
            expected(1:end-1, 3:end) = plain(2:end, 1:end-2);
            testCase.verifyEqual(isnan(moved), isnan(expected));
            testCase.verifyEqual(moved(~isnan(moved)), expected(~isnan(expected)), 'RelTol', single(1e-6));
        end

        function testSwitchWithinWindow(testCase)
            S2 = testCase.makeData();
            n = testCase.totalLines;
            x = zeros(n, 1); y = zeros(n, 1);
            x(53:end) = 1; y(57:end) = -2;      % switches in the middle of cycles
            off = testCase.offsets(x, y);
            map = [0 1; 1 0];
            [IM, F] = slap2Frames(S2, [1 2], [50 55 60], 3, CompensateOnline = true, Offsets = off, RefOffset = [0 0]);
            [IMref, Fref] = testCase.reference(S2, [1 2], [50 55 60], 3, 'exp', [x y], map);
            testCase.verifyEqual(isnan(IM), isnan(IMref));
            testCase.verifyEqual(double(IM(~isnan(IM))), IMref(~isnan(IMref)), 'RelTol', 1e-6);
            testCase.verifyEqual(double(F), Fref, 'RelTol', 1e-6);
        end

        function testCycleKernel(testCase)
            S2 = testCase.makeData();
            [IM, F] = slap2Frames(S2, 1, [1 4 20], 1, Kernel = 'cycle');
            [IMref, Fref] = testCase.reference(S2, 1, [1 4 20], 1, 'cycle', zeros(testCase.totalLines, 2), [0 1; 1 0]);
            testCase.verifyEqual(isnan(IM), isnan(IMref));
            testCase.verifyEqual(double(IM(~isnan(IM))), IMref(~isnan(IMref)), 'RelTol', 1e-6);
            testCase.verifyEqual(double(F), Fref, 'RelTol', 1e-6);
        end

        function testRowsAndColsWindow(testCase)
            S2 = testCase.makeData();
            full = slap2Frames(S2, [1 2], [20 40], 4);
            part = slap2Frames(S2, [1 2], [20 40], 4, Rows = 3:8, Cols = 2:10);
            expected = full(3:8, 2:10, :, :);
            testCase.verifyEqual(isnan(part), isnan(expected));
            testCase.verifyEqual(part(~isnan(part)), expected(~isnan(expected)), 'RelTol', single(1e-6));
        end

        function testInvalidSamplesAreSkipped(testCase)
            S2 = testCase.makeData(true);
            [IM, F] = slap2Frames(S2, [1 2], [30 31], 4);
            [IMref, Fref] = testCase.reference(S2, [1 2], [30 31], 4, 'exp', zeros(testCase.totalLines, 2), [0 1; 1 0]);
            testCase.verifyEqual(double(IM(~isnan(IM))), IMref(~isnan(IMref)), 'RelTol', 1e-6);
            testCase.verifyEqual(double(F), Fref, 'RelTol', 1e-6);
        end

        function testNonIntegerLinesError(testCase)
            S2 = testCase.makeData();
            testCase.verifyError(@() slap2Frames(S2, 1, 20.5, 4), 'slap2Frames:nonIntegerLines');
        end

        function testGetImagesSemantics(testCase)
            % Boundary 'exclusive' and MaskInvalid false give getImages' frames:
            % lines exactly dtRead away are left out and invalid samples are summed
            S2 = testCase.makeData(true);
            times = [20 30 31 70];
            [Num, Den, pix] = slap2FrameSums(S2, [2 1], times, 4, Boundary = 'exclusive', MaskInvalid = false);
            [IMref, Fref] = testCase.reference(S2, [2 1], times, 4, 'exp', zeros(testCase.totalLines, 2), [0 1; 1 0], true, false);
            [IM, F] = testCase.toImage(Num, Den, pix, numel(times));
            testCase.verifySize(Den, [numel(pix) numel(times)]);
            testCase.verifyEqual(isnan(IM), isnan(IMref));
            testCase.verifyEqual(IM(~isnan(IM)), IMref(~isnan(IMref)), 'RelTol', 1e-9);
            testCase.verifyEqual(F, Fref, 'RelTol', 1e-9);
        end

        function testPixelsSubset(testCase)
            S2 = testCase.makeData();
            n = testCase.totalLines;
            x = zeros(n, 1); y = zeros(n, 1); x(53:end) = 1; y(57:end) = -2;
            args = {'CompensateOnline', true, 'Offsets', testCase.offsets(x, y), 'RefOffset', [0 0]};
            [NumAll, DenAll, pixAll] = slap2FrameSums(S2, [1 2], [50 55 60], 3, args{:});
            want = [5 17 18 40 41 66 101 119 120]';
            [Num, Den, pix] = slap2FrameSums(S2, [1 2], [50 55 60], 3, args{:}, 'Pixels', want);
            [tf, loc] = ismember(pix, pixAll);
            testCase.verifyTrue(all(tf));
            testCase.verifyEqual(pix, intersect(want, pixAll));
            testCase.verifyEqual(Num, NumAll(loc, :, :), 'RelTol', 1e-12);
            testCase.verifyEqual(Den, DenAll(loc, :, :), 'RelTol', 1e-12);
        end

        function testFileReaderMatchesGetLineData(testCase)
            import matlab.unittest.fixtures.TemporaryFolderFixture
            folder = testCase.applyFixture(TemporaryFolderFixture).Folder;
            S2 = testCase.makeDiskData(folder, true);
            n = testCase.totalLines;
            x = zeros(n, 1); y = zeros(n, 1); x(53:end) = 1; y(57:end) = -2; x(80:end) = -1;
            off = testCase.offsets(x, y);
            times = [20 50 55 60 72 75 115];   % blocks span both files
            for exclusive = [false true]
                bnd = 'inclusive'; if exclusive, bnd = 'exclusive'; end
                args = {[2 1], times, 3, 'Boundary', bnd, 'MaskInvalid', ~exclusive, ...
                    'CompensateOnline', true, 'Offsets', off, 'RefOffset', [0 0]};
                [Nf, Df, pf] = slap2FrameSums(S2, args{:}, 'Reader', 'file');
                [Ng, Dg, pg] = slap2FrameSums(S2, args{:}, 'Reader', 'getLineData');
                testCase.verifyEqual(pf, pg);
                testCase.verifyEqual(Nf, Ng);
                testCase.verifyEqual(Df, Dg);
                [IM, F] = testCase.toImage(Nf, Df, pf, numel(times));
                [IMref, Fref] = testCase.reference(S2, [2 1], times, 3, 'exp', [x y], [0 1; 1 0], exclusive, ~exclusive);
                testCase.verifyEqual(isnan(IM), isnan(IMref));
                testCase.verifyEqual(IM(~isnan(IM)), IMref(~isnan(IMref)), 'RelTol', 1e-9);
                testCase.verifyEqual(F, Fref, 'RelTol', 1e-9);
            end
            % the default reader is the file reader, also right after another recording was read
            other = testCase.makeDiskData(folder, false, 'other');
            slap2Frames(other, 1, 40, 4);
            [IM1, F1] = slap2Frames(S2, 1, 40, 4);
            [IMref, Fref] = testCase.reference(S2, 1, 40, 4, 'exp', zeros(n, 2), [0 1; 1 0]);
            testCase.verifyEqual(isnan(IM1), isnan(IMref));
            testCase.verifyEqual(double(IM1(~isnan(IM1))), IMref(~isnan(IMref)), 'RelTol', 1e-6);
            testCase.verifyEqual(double(F1), Fref, 'RelTol', 1e-6);
        end
    end

    methods (Access = private)
        function n = totalLines(testCase)
            n = testCase.P * testCase.nCycles;
        end

        function S2 = makeData(testCase, withInvalid)
            if nargin < 2, withInvalid = false; end
            nPx = testCase.NR * testCase.NC;
            % each line samples a band of raster pixels, some twice per cycle, plus one block superpixel
            ids = cell(1, testCase.P);
            for L = 1:testCase.P
                k = (L-1)*20 + (1:25);
                k = mod(k - 1, nPx) + 1;
                ids{L} = uint32([k(:) - 1; nPx + L - 1]);   % 0-based; the last id is not raster
            end
            M = struct();
            M.header = struct('linesPerCycle', testCase.P, 'dmdPixelsPerRow', testCase.NR, ...
                'dmdPixelsPerColumn', testCase.NC);
            M.totalNumLines = testCase.totalLines;
            M.numCycles = testCase.nCycles;
            M.lineFastZIdxs = ones(1, testCase.P);
            M.lineSuperPixelIDs = ids;
            M.getLineData = @(li, ci, ch) TestSlap2Frames.lineData(ids, li, ci, ch, withInvalid);
            S2 = struct('hMultiDataFiles', M);
        end

        function S2 = makeDiskData(testCase, folder, withInvalid, tag)
            % makeData's recording, also written to two data files: a file
            % header, then whole cycles of lines (line header, channel 1 data,
            % channel 2 data) with padding at the end of each cycle
            if nargin < 4, tag = 'rec'; end
            S2 = testCase.makeData(withInvalid);
            M = S2.hMultiDataFiles;
            P = testCase.P; nCh = 2; H = 10; hw = 4; pad = 2;
            nSp = cellfun(@numel, M.lineSuperPixelIDs(:));
            lineWords = hw + nCh * nSp;
            lineHeaderIdxs = H + 1 + [0; cumsum(lineWords(1:end-1))];
            cycleWords = sum(lineWords) + pad;
            M.filename = fullfile(folder, [tag '.dat']);
            M.header.numChannels = nCh;
            M.header.channelsInterleave = false;
            M.header.bytesPerCycle = 2 * cycleWords;
            M.lineHeaderIdxs = lineHeaderIdxs;
            M.lineDataStartIdxs = lineHeaderIdxs + hw;
            M.lineDataNumElements = nCh * nSp;
            fileCycles = {1:12, 13:testCase.nCycles};
            files = struct('datFileName', {}, 'header', {}, 'lineHeaderIdxs', {}, 'numCycles', {});
            for k = 1:numel(fileCycles)
                words = zeros(H + cycleWords * numel(fileCycles{k}), 1, 'int16');
                words(1:H) = int16(-(1:H));
                for j = 1:numel(fileCycles{k})
                    c = fileCycles{k}(j);
                    for L = 1:P
                        d = TestSlap2Frames.lineData(M.lineSuperPixelIDs, L, c, 1:nCh, withInvalid);
                        at = (j - 1) * cycleWords + lineHeaderIdxs(L);
                        words(at + (0:hw-1)) = int16(-100 - L);
                        words(at + hw + (0:numel(d{1}) - 1)) = d{1}(:);
                    end
                end
                fn = fullfile(folder, sprintf('%s-CYCLE-%02d.dat', tag, k));
                fid = fopen(fn, 'w'); fwrite(fid, words, 'int16'); fclose(fid);
                files(k) = struct('datFileName', fn, 'header', M.header, ...
                    'lineHeaderIdxs', lineHeaderIdxs, 'numCycles', numel(fileCycles{k}));
            end
            M.hDataFiles = files;
            S2 = struct('hMultiDataFiles', M);
        end

        function [IM, F] = toImage(testCase, Num, Den, pix, nF)
            % images (rows x cols x channels x frames) from slap2FrameSums' sums
            NR = testCase.NR; NC = testCase.NC; nCh = size(Num, 3);
            t = floor((pix - 1) / NR) + 1 + NC * mod(pix - 1, NR);
            IM = nan(NC * NR, nF, nCh); F = zeros(NC * NR, nF);
            IM(t, :, :) = Num ./ Den(:, :, 1);
            F(t, :) = Den(:, :, 1);
            IM = permute(reshape(IM, NC, NR, nF, nCh), [1 2 4 3]);
            F = reshape(F, NC, NR, nF);
        end

        function [IM, F] = reference(testCase, S2, channels, times, dt, kernel, xy, map, exclusive, mask)
            % exclusive: leave out lines exactly dtRead away (getImages); mask:
            % skip invalid samples (getImage)
            if nargin < 9, exclusive = false; end
            if nargin < 10, mask = true; end
            M = S2.hMultiDataFiles; P = testCase.P; NR = testCase.NR; NC = testCase.NC; nPx = NR * NC;
            n = testCase.totalLines;
            IM = nan(NC, NR, numel(channels), numel(times)); F = zeros(NC, NR, numel(times));
            for f = 1:numel(times)
                if strcmp(kernel, 'exp')
                    dtRead = max(3*dt, P);
                    lns = max(1, times(f) - dtRead):min(n, times(f) + dtRead);
                    if exclusive, lns = lns(abs(times(f) - lns) < dtRead); end
                    w = exp(-abs((times(f) - lns) / dt));
                else
                    lns = ((times(f)-1)*P + 1):(times(f)*P); w = ones(size(lns));
                end
                for k = 1:numel(channels)
                    num = zeros(NC, NR); den = zeros(NC, NR);
                    for j = 1:numel(lns)
                        l = lns(j); L = mod(l-1, P) + 1; c = floor((l-1)/P) + 1;
                        d = M.getLineData(L, c, channels(k)); d = double(d{1}(:, 1));
                        ids = double(M.lineSuperPixelIDs{L}(:)) + 1;
                        shift = map * xy(l, :)';
                        for s = 1:numel(ids)
                            if ids(s) > nPx || (mask && d(s) <= double(intmin('int16'))), continue, end
                            r = floor((ids(s)-1)/NR) + 1 + shift(1); cc = mod(ids(s)-1, NR) + 1 + shift(2);
                            if r < 1 || r > NC || cc < 1 || cc > NR, continue, end
                            num(r, cc) = num(r, cc) + w(j) * d(s);
                            den(r, cc) = den(r, cc) + w(j);
                        end
                    end
                    IM(:, :, k, f) = num ./ den;
                    if k == 1, F(:, :, f) = den; end
                end
            end
        end
    end

    methods (Static)
        function off = offsets(x, y)
            off = struct('x', int16(x), 'y', int16(y), 'z', zeros(numel(x), 1, 'single'));
        end

        function out = lineData(ids, li, ci, ch, withInvalid)
            % channels min(ch):max(ch) in file order, as slap2.util.DataFile.getLineData does
            chs = min(ch):max(ch);
            out = cell(numel(li), 1);
            for i = 1:numel(li)
                k = double(ids{li(i)}(:));
                v = zeros(numel(k), numel(chs));
                for j = 1:numel(chs)
                    v(:, j) = 1000 * chs(j) + 3 * k + 7 * li(i) + 11 * ci(i) + mod(k .* ci(i), 5);
                end
                if withInvalid && mod(ci(i), 3) == 0
                    v(2, :) = double(intmin('int16'));
                end
                out{i} = int16(v);
            end
        end
    end
end
