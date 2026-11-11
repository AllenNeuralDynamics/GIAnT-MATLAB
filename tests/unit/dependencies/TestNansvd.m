classdef TestNansvd < matlab.unittest.TestCase
    %TESTNANSVD Unit tests for nansvd.
    %   nansvd drops pixels (rows) and frames (columns) with too many NaNs,
    %   subtracts a moving-mean background, and fits a rank-numComps SVD,
    %   filling the remaining NaNs iteratively. svds returns
    %   min(numComps, min(size(D))) components, so a user ROI with fewer
    %   valid pixels, or a trial with fewer frames, than numComps gets
    %   fewer; nansvd pads the outputs to numComps components (NaN in U
    %   and V, zero singular values in S).
    %   Run with: runtests('tests/unit/dependencies/TestNansvd.m')

    properties (Constant)
        NumComps = 3
        MaxIter = 10
        NanFrac = 0.3
    end

    properties (TestParameter)
        nValidPx = struct('onePixel', 1, 'twoPixels', 2)
    end

    methods (TestClassSetup)
        function addPaths(testCase)
            import matlab.unittest.fixtures.PathFixture
            testsDir = fileparts(mfilename('fullpath'));
            repoRoot = fileparts(fileparts(fileparts(testsDir)));
            testCase.applyFixture(PathFixture(fullfile(repoRoot, 'dependencies')));
        end
    end

    methods (TestMethodSetup)
        function seedRng(testCase)
            testCase.addTeardown(@rng, rng);
            rng(0, 'twister');
        end
    end

    methods (Test)
        function testFewPixelsPadsToRequestedSize(testCase, nValidPx)
            % 4525 valid frames, as in the user ROI that failed on the
            % 780125 Neuron1 dataset, plus unscanned pixels and frames.
            [D, validPx, validFr] = makeRoi(nValidPx, 4525, 4, 7);
            k = testCase.NumComps;

            [U, S, V, BG] = nansvd(D, k, testCase.MaxIter, testCase.NanFrac);

            testCase.verifySize(U, [size(D,1) k]);
            testCase.verifySize(S, [k k]);
            testCase.verifySize(V, [size(D,2) k]);
            testCase.verifySize(BG, size(D));

            % components svds did not return are padded
            testCase.verifyTrue(all(isnan(U(:, nValidPx+1:end)), 'all'));
            testCase.verifyTrue(all(isnan(V(:, nValidPx+1:end)), 'all'));
            Spad = S;
            Spad(1:nValidPx, 1:nValidPx) = 0;
            testCase.verifyEqual(Spad, zeros(k));

            % the returned components cover the valid pixels and frames only
            testCase.verifyTrue(all(isfinite(U(validPx, 1:nValidPx)), 'all'));
            testCase.verifyTrue(all(isfinite(V(validFr, 1:nValidPx)), 'all'));
            testCase.verifyTrue(all(isnan(U(~validPx, :)), 'all'));
            testCase.verifyTrue(all(isnan(V(~validFr, :)), 'all'));
        end

        function testPaddingKeepsReturnedComponents(testCase, nValidPx)
            % Asking for more components than svds can return gives the
            % same components as asking for exactly as many as it returns.
            D = makeRoi(nValidPx, 600, 4, 7);
            D(rand(size(D)) < 0.05) = NaN; % gaps to fill iteratively

            rng(1, 'twister');
            [Uk, Sk, Vk, BGk] = nansvd(D, nValidPx, testCase.MaxIter, testCase.NanFrac);
            rng(1, 'twister');
            [U, S, V, BG] = nansvd(D, testCase.NumComps, testCase.MaxIter, testCase.NanFrac);

            testCase.verifyEqual(U(:, 1:nValidPx), Uk);
            testCase.verifyEqual(S(1:nValidPx, 1:nValidPx), Sk);
            testCase.verifyEqual(V(:, 1:nValidPx), Vk);
            testCase.verifyEqual(BG, BGk);
        end

        function testShortTrialPadsToRequestedSize(testCase)
            % A two-frame trial has fewer timepoints than components.
            D = 100 + 5*randn(40, 2);
            k = testCase.NumComps;

            [U, S, V, BG] = nansvd(D, k, testCase.MaxIter, testCase.NanFrac);

            testCase.verifySize(U, [40 k]);
            testCase.verifySize(S, [k k]);
            testCase.verifySize(V, [2 k]);
            testCase.verifySize(BG, [40 2]);
            testCase.verifyTrue(all(isnan(U(:, 3:end)), 'all'));
            testCase.verifyTrue(all(isnan(V(:, 3:end)), 'all'));
        end

        function testAllComponentsKeepSingularValues(testCase)
            % When svds returns every requested component, nothing is
            % padded and no singular value is overwritten.
            [D, validPx, validFr] = makeRoi(40, 600, 4, 7);
            D(rand(size(D)) < 0.05) = NaN;
            k = testCase.NumComps;

            [U, S, V] = nansvd(D, k, testCase.MaxIter, testCase.NanFrac);

            testCase.verifySize(U, [size(D,1) k]);
            testCase.verifySize(S, [k k]);
            testCase.verifySize(V, [size(D,2) k]);
            s = diag(S);
            testCase.verifyEqual(S, diag(s));
            testCase.verifyGreaterThan(s, 0);
            testCase.verifyTrue(issorted(s, 'descend'));
            testCase.verifyTrue(all(isfinite(U(validPx, :)), 'all'));
            testCase.verifyTrue(all(isfinite(V(validFr, :)), 'all'));
        end
    end
end

function [D, validPx, validFr] = makeRoi(nPx, nFrames, nUnscannedPx, nUnscannedFr)
% nPx pixels sharing a transient-like signal plus noise, with nUnscannedPx
% all-NaN pixels and nUnscannedFr all-NaN frames mixed in, which nansvd's
% NaN-fraction selection drops.
t = 1:nFrames;
sig = 1 + 0.5*exp(-mod(t, 300)/40);
D = 100*(rand(nPx, 1) + 0.5)*sig + 5*randn(nPx, nFrames);
D = [D; nan(nUnscannedPx, nFrames)];
D(:, randperm(nFrames, nUnscannedFr)) = NaN;
order = randperm(size(D, 1));
D = D(order, :);
validPx = order(:) <= nPx;
validFr = ~all(isnan(D), 1)';
end
