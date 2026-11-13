classdef TestGetLineOnlineOffsets < matlab.unittest.TestCase
    %TESTGETLINEONLINEOFFSETS Unit tests for getLineOnlineOffsets.
    %   Builds line headers in an int16 word array laid out like a SLAP2
    %   .dat file (bytes 16-19 flags, 20-21 x, 22-23 y, 32-35 z) and checks
    %   that every line decodes, across several files.
    %   Run with: runtests('tests/unit/dependencies/alignment/slap2/TestGetLineOnlineOffsets.m')

    methods (TestClassSetup)
        function addPaths(testCase)
            import matlab.unittest.fixtures.PathFixture
            testsDir = fileparts(mfilename('fullpath'));
            repoRoot = fileparts(fileparts(fileparts(fileparts(fileparts(testsDir)))));
            testCase.applyFixture(PathFixture(fullfile(repoRoot, 'dependencies', 'alignment', 'slap2')));
        end
    end

    methods (Test)
        function testDecodesEveryLine(testCase)
            P = 4; nC = 3; n = P * nC;
            x = int16((1:n)' - 5); y = int16(-(1:n)'); z = single((1:n)' / 10);
            pend = mod((1:n)', 3) == 0; pin = mod((1:n)', 4) < 2;
            D = TestGetLineOnlineOffsets.makeFile(P, x, y, z, pend, pin, 40);
            off = getLineOnlineOffsets(D);
            testCase.verifyEqual(off.x, x);
            testCase.verifyEqual(off.y, y);
            testCase.verifyEqual(off.z, z);
            testCase.verifyEqual(off.syncPending, pend);
            testCase.verifyEqual(off.dmdSyncPin, pin);
            testCase.verifyEqual(off.linesPerCycle, P);
            testCase.verifyEqual(off.numLines, n);
        end

        function testConcatenatesFilesInOrder(testCase)
            P = 3;
            x1 = int16([0 0 0 1 1 1]'); x2 = int16([1 2 2 2 2 2]');
            y1 = zeros(6, 1, 'int16'); y2 = int16([0 0 0 0 -1 -1]');
            f = false(6, 1);
            D1 = TestGetLineOnlineOffsets.makeFile(P, x1, y1, zeros(6, 1, 'single'), f, f, 40);
            D2 = TestGetLineOnlineOffsets.makeFile(P, x2, y2, zeros(6, 1, 'single'), f, f, 40);
            off = getLineOnlineOffsets(struct('hMultiDataFiles', struct('hDataFiles', [D1 D2])));
            testCase.verifyEqual(off.x, [x1; x2]);
            testCase.verifyEqual(off.xySwitchLines, [4; 8; 11]);
        end

        function testDmdYTiming(testCase)
            % a new y reaches the DMD 2 lines before the header records it, or
            % 6 lines when x changes too, but not before the previous switch
            P = 10; n = 70;
            x = zeros(n, 1, 'int16'); y = zeros(n, 1, 'int16');
            y(20:end) = 1;                 % y only at line 20
            x(40:end) = 1; y(40:end) = 2;  % x and y at line 40
            x(50:end) = 2;                 % x only at line 50
            y(51:end) = 3;                 % y only at 51, one line after the previous switch
            f = false(n, 1);
            D = TestGetLineOnlineOffsets.makeFile(P, x, y, zeros(n, 1, 'single'), f, f, 40);
            off = getLineOnlineOffsets(D);
            expected = y;
            expected(18:19) = 1;
            expected(34:39) = 2;
            expected(50) = 3;
            testCase.verifyEqual(off.y, y);        % the header values are unchanged
            testCase.verifyEqual(off.yEff, expected);
            testCase.verifyEqual(off.yLead, [2 6]);
            off3 = getLineOnlineOffsets(D, YLead = [3 0]);
            expected = y; expected(17:19) = 1; expected(50) = 3;
            testCase.verifyEqual(off3.yEff, expected);
        end

        function testHeadersWithoutZ(testCase)
            P = 2; n = 4;
            D = TestGetLineOnlineOffsets.makeFile(P, int16((1:n)'), int16((1:n)'), single(ones(n, 1)), false(n, 1), false(n, 1), 32);
            off = getLineOnlineOffsets(D);
            testCase.verifyEqual(off.z, zeros(n, 1, 'single'));
        end
    end

    methods (Static, Access = private)
        function D = makeFile(P, x, y, z, pend, pin, headerBytes)
            nC = numel(x) / P;
            wordsPerLine = headerBytes / 2 + 8;            % header plus some line data
            wordsPerCycle = P * wordsPerLine;
            data = zeros(wordsPerCycle * nC + 10, 1, 'int16');
            lineHeaderIdxs = 5 + (0:P-1) * wordsPerLine;   % 1-based word index of each header
            for c = 1:nC
                for L = 1:P
                    i = (c-1)*P + L;
                    base = lineHeaderIdxs(L) + (c-1) * wordsPerCycle;
                    flags = uint32(bitshift(uint32(pend(i)), 5) + bitshift(uint32(pin(i)), 6));
                    data(base + (8:9)) = typecast(flags, 'int16');
                    data(base + 10) = x(i); data(base + 11) = y(i);
                    if headerBytes > 32
                        data(base + (16:17)) = typecast(single(z(i)), 'int16');
                    end
                end
            end
            D = struct('rawData', struct('Format', 'int16', 'Data', data), ...
                'lineHeaderIdxs', lineHeaderIdxs, 'numCycles', nC, ...
                'header', struct('bytesPerCycle', 2 * wordsPerCycle, 'lineHeaderSizeBytes', headerBytes));
        end
    end
end
