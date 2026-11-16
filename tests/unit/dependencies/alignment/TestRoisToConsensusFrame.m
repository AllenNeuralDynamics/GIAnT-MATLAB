classdef TestRoisToConsensusFrame < matlab.unittest.TestCase
    %TESTROISTOCONSENSUSFRAME Unit tests for roisToConsensusFrame.
    %   ROIs drawn on an annotation trial move by -motOutput of that trial
    %   onto the consensus frame, or through its full map when SILo's warp is
    %   affine or patchwise. When the trial is unknown or unaligned they stay
    %   put, but their masks still take the consensus frame size.
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

        function testShiftWarpMatchesOffset(testCase)
            % a rigid warp takes the same path as the offsets alone
            entry = testCase.entry('E1T2DMD1_REGISTERED_DOWNSAMPLED-20Hz', false(10));
            warp = testCase.affineWarp([1 0 2; 0 1 -1]);

            out = roisToConsensusFrame(entry, testCase.FnRegDS, testCase.MotOutput, testCase.TargetSz, warp);

            testCase.verifyEqual(find(out.roiData{1}.mask), sub2ind(testCase.TargetSz, 3, 6));
            testCase.verifyEqual(out.roiData{1}.Center, [6 3]);
        end

        function testAffineWarpMapsThroughFullMap(testCase)
            % rows stretched 10% about the origin: common (r, c) is read from the
            % trial at (1.1 r, c). Trial pixel (33, 10) lands on common (30, 10);
            % a shift by the displacement at the centre (2.05 rows) would put it
            % on row 31
            sz = [40 40];
            mask = false(sz); mask(33, 10) = true;
            roi = struct('Label', 'ROI', 'mask', mask, 'Center', [10 33]);
            entry = struct('fn', 'E1T2DMD1_REGISTERED_DOWNSAMPLED-20Hz', 'roiData', {{roi}});
            warp = testCase.affineWarp([1.1 0 0; 0 1 0]);
            motOutput = [0 2.05; 0 0];

            out = roisToConsensusFrame(entry, testCase.FnRegDS, motOutput, sz, warp);

            testCase.verifyEqual(find(out.roiData{1}.mask), sub2ind(sz, 30, 10));
            testCase.verifyEqual(out.roiData{1}.Center, [10 30], 'AbsTol', 1e-9);
        end

        function testPatchwiseWarpMapsThroughField(testCase)
            % the same stretch as a patchwise field (which reproduces an affine
            % map exactly) gives the same ROI
            sz = [40 40];
            mask = false(sz); mask(33, 10) = true;
            roi = struct('Label', 'ROI', 'mask', mask, 'Center', [10 33]);
            entry = struct('fn', 'E1T2DMD1_REGISTERED_DOWNSAMPLED-20Hz', 'roiData', {{roi}});
            A = [1.1 0 0; 0 1 0];
            [nc, nr] = meshgrid(5:10:35, 5:10:35); nodes = [nr(:) nc(:)];
            d = A(:, 1:2)*nodes' + A(:, 3) - nodes';
            warp = struct('affine', cat(3, [1 0 0; 0 1 0], A), 'nodes', nodes, 'sigma', 10, ...
                'disp', cat(3, zeros(size(d)), d), 'patch', []);

            out = roisToConsensusFrame(entry, testCase.FnRegDS, [0 2.05; 0 0], sz, warp);

            testCase.verifyEqual(find(out.roiData{1}.mask), sub2ind(sz, 30, 10));
            testCase.verifyEqual(out.roiData{1}.Center, [10 30], 'AbsTol', 1e-6);
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

        function warp = affineWarp(~, A2)
            % trial 1 identity, trial 2 the affine map A2
            warp = struct('affine', cat(3, [1 0 0; 0 1 0], A2), 'nodes', zeros(0, 2), ...
                'sigma', nan, 'disp', zeros(2, 0, 2), 'patch', []);
        end
    end
end
