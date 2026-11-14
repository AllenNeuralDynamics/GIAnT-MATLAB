classdef TestRoisToConsensusFrame < matlab.unittest.TestCase
    %TESTROISTOCONSENSUSFRAME Unit tests for roisToConsensusFrame.
    %   ROIs drawn on an annotation trial move by -motOutput of that trial
    %   onto the consensus frame. When the trial is unknown or unaligned
    %   they stay put, but their masks still take the consensus frame size.
    %   Run with: runtests('tests/unit/dependencies/alignment/TestRoisToConsensusFrame.m')

    properties
        FnRegDS = {'E1T1DMD1_REGISTERED_DOWNSAMPLED-20Hz.h5', 'E1T2DMD1_REGISTERED_DOWNSAMPLED-20Hz.h5'}
        MotOutput = [0 2; 0 -1]   % trial 2: aligned(r,c) = trial(r+2, c-1)
        TargetSz = [12 12]
    end

    methods (TestClassSetup)
        function addPaths(testCase)
            import matlab.unittest.fixtures.PathFixture
            testsDir = fileparts(mfilename('fullpath'));
            repoRoot = fileparts(fileparts(fileparts(fileparts(testsDir))));
            testCase.applyFixture(PathFixture(fullfile(repoRoot, 'dependencies', 'alignment')));
        end
    end

    methods (Test)
        function testShiftsByAnnotationTrialOffset(testCase)
            entry = testCase.entry('E1T2DMD1_REGISTERED_DOWNSAMPLED-20Hz', false(10));

            out = roisToConsensusFrame(entry, testCase.FnRegDS, testCase.MotOutput, testCase.TargetSz);

            testCase.verifyEqual(find(out.roiData{1}.mask), sub2ind(testCase.TargetSz, 3, 6));
            testCase.verifyEqual(out.roiData{1}.Center, [6 3]); % [x y]: col 5 -> 6, row 5 -> 3
        end

        function testUnknownTrialKeepsPositionOnConsensusCanvas(testCase)
            entry = testCase.entry('not_a_trial', false(15));

            out = testCase.verifyWarning(@() roisToConsensusFrame(entry, testCase.FnRegDS, ...
                testCase.MotOutput, testCase.TargetSz), 'roisToConsensusFrame:UnknownTrial');

            testCase.verifySize(out.roiData{1}.mask, testCase.TargetSz);
            testCase.verifyEqual(find(out.roiData{1}.mask), sub2ind(testCase.TargetSz, 5, 5));
            testCase.verifyEqual(out.roiData{1}.Center, [5 5]);
        end

        function testUnalignedTrialKeepsPositionOnConsensusCanvas(testCase)
            entry = testCase.entry('E1T2DMD1_REGISTERED_DOWNSAMPLED-20Hz', false(8));
            motOutput = testCase.MotOutput;
            motOutput(:, 2) = NaN;

            out = testCase.verifyWarning(@() roisToConsensusFrame(entry, testCase.FnRegDS, ...
                motOutput, testCase.TargetSz), 'roisToConsensusFrame:TrialNotAligned');

            testCase.verifySize(out.roiData{1}.mask, testCase.TargetSz);
            testCase.verifyEqual(find(out.roiData{1}.mask), sub2ind(testCase.TargetSz, 5, 5));
        end

        function testEntryWithoutRoisIsUnchanged(testCase)
            entry = struct('fn', 'not_a_trial', 'roiData', {{}});

            out = roisToConsensusFrame(entry, testCase.FnRegDS, testCase.MotOutput, testCase.TargetSz);

            testCase.verifyEqual(out, entry);
        end
    end

    methods
        function entry = entry(~, fn, mask)
            % one ROI covering pixel (5,5) of a mask of the given size
            mask(5, 5) = true;
            roi = struct('Label', 'ROI', 'mask', mask, 'Center', [5 5]);
            entry = struct('fn', fn, 'roiData', {{roi}});
        end
    end
end
