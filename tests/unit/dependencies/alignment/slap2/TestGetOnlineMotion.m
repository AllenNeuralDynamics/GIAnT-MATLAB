classdef TestGetOnlineMotion < matlab.unittest.TestCase
    %TESTGETONLINEMOTION Unit tests for getOnlineMotion.
    %   Uses a struct in place of the SLAP2 data file. Its getLineHeader
    %   encodes the requested (line, cycle) in the returned offsets:
    %   xOffset_pix = line, yOffset_pix = cycle, zOffset_um = 1000*cycle + line.
    %   Run with: runtests('tests/unit/dependencies/alignment/slap2/TestGetOnlineMotion.m')

    methods (TestClassSetup)
        function addPaths(testCase)
            import matlab.unittest.fixtures.PathFixture
            testsDir = fileparts(mfilename('fullpath'));
            repoRoot = fileparts(fileparts(fileparts(fileparts(fileparts(testsDir)))));
            testCase.applyFixture(PathFixture(fullfile(repoRoot, 'dependencies', 'alignment', 'slap2')));
        end
    end

    methods (Test)
        function testDefaultLineIdsOnePerCycle(testCase)
            hDataFile = testCase.makeDataFile(4, 3);

            [x, y, z] = getOnlineMotion(hDataFile);

            % Defaults to the first line of every cycle.
            testCase.verifyEqual(x, [1; 1; 1]);
            testCase.verifyEqual(y, [1; 2; 3]);
            testCase.verifyEqual(z, [1001; 2001; 3001]);
        end

        function testLineIdsMapToLineAndCycle(testCase)
            hDataFile = testCase.makeDataFile(4, 3);

            [x, y, z] = getOnlineMotion(hDataFile, [1 4 5 12]);

            testCase.verifyEqual(x, [1; 4; 1; 4]);
            testCase.verifyEqual(y, [1; 1; 2; 3]);
            testCase.verifyEqual(z, [1001; 1004; 2001; 3004]);
        end

        function testColumnLineIds(testCase)
            hDataFile = testCase.makeDataFile(4, 3);

            [x, y, z] = getOnlineMotion(hDataFile, [2; 7]);

            testCase.verifyEqual(x, [2; 3]);
            testCase.verifyEqual(y, [1; 2]);
            testCase.verifyEqual(z, [1002; 2003]);
        end

        function testLinesPerCycleFromHeaderCount(testCase)
            % Only the number of line headers matters, not their values.
            hDataFile = testCase.makeDataFile(3, 2);
            hDataFile.lineHeaderIdxs = [10 20 30];

            [x, y] = getOnlineMotion(hDataFile);

            testCase.verifyEqual(x, [1; 1]);
            testCase.verifyEqual(y, [1; 2]);
        end

        function testEmptyLineIds(testCase)
            hDataFile = testCase.makeDataFile(4, 3);

            [x, y, z] = getOnlineMotion(hDataFile, []);

            testCase.verifySize(x, [0 1]);
            testCase.verifySize(y, [0 1]);
            testCase.verifySize(z, [0 1]);
        end
    end

    methods (Static, Access = private)
        function hDataFile = makeDataFile(linesPerCycle, numCycles)
            hDataFile = struct( ...
                'lineHeaderIdxs', 1:linesPerCycle, ...
                'numCycles', numCycles, ...
                'getLineHeader', @(line, cycle) struct('xOffset_pix', line, ...
                    'yOffset_pix', cycle, 'zOffset_um', 1000*cycle + line));
        end
    end
end
