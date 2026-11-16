classdef TestFrameOnlineOffsets < matlab.unittest.TestCase
    %TESTFRAMEONLINEOFFSETS Unit tests for frameOnlineOffsets.
    %   Run with: runtests('tests/unit/dependencies/alignment/slap2/TestFrameOnlineOffsets.m')

    methods (TestClassSetup)
        function addPaths(testCase)
            import matlab.unittest.fixtures.PathFixture
            testsDir = fileparts(mfilename('fullpath'));
            repoRoot = fileparts(fileparts(fileparts(fileparts(fileparts(testsDir)))));
            testCase.applyFixture(PathFixture(fullfile(repoRoot, 'dependencies', 'alignment', 'slap2')));
        end
    end

    methods (Test)
        function testKernelWeightedRoundedMean(testCase)
            n = 400; P = 10; dt = 20;
            x = zeros(n, 1); y = zeros(n, 1);
            x(201:end) = 3; y(201:end) = -2; yEff = y; yEff(195:200) = -2;
            off = struct('x', int16(x), 'y', int16(y), 'yEff', int16(yEff), 'linesPerCycle', P);
            times = [100 196 205 300];
            o = frameOnlineOffsets(off, times, dt);
            expected = zeros(numel(times), 2);
            for f = 1:numel(times)
                lns = max(1, times(f) - 60):min(n, times(f) + 60);
                w = exp(-abs(times(f) - lns(:)) / dt);
                expected(f, :) = round([w' * x(lns), w' * yEff(lns)] / sum(w));
            end
            testCase.verifyEqual(o, expected);
            testCase.verifyEqual(o(1, :), [0 0]);
            testCase.verifyEqual(o(end, :), [3 -2]);
        end

        function testFallsBackToHeaderY(testCase)
            off = struct('x', int16([0 0 1 1]'), 'y', int16([0 0 1 1]'), 'linesPerCycle', 2);
            o = frameOnlineOffsets(off, 4, 1);
            testCase.verifyEqual(o, [1 1]);
        end
    end
end
