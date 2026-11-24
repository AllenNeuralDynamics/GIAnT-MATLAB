classdef TestDebiasVariance < matlab.unittest.TestCase
    %TESTDEBIASVARIANCE params.debiasVariance, the variance term of
    %extraction's final (debias) refit of the activity (issue #56): 'fixed'
    %by default, also when the field is missing; 'model' keeps the earlier
    %fit. Only the activity depends on it: footprints, F0 and the
    %least-squares activity are identical. The amplitude itself is tested
    %in TestSubproblemRefitTest/testDebiasRecoversFoundEvents.
    %   Run with: runtests('tests/unit/dependencies/extraction/TestDebiasVariance.m')

    properties
        Inputs
    end

    methods (TestClassSetup)
        function setup(testCase)
            import matlab.unittest.fixtures.PathFixture
            import matlab.unittest.fixtures.TemporaryFolderFixture
            testsDir = fileparts(mfilename('fullpath'));
            repoRoot = fileparts(fileparts(fileparts(fileparts(testsDir))));
            testCase.applyFixture(PathFixture(fullfile(repoRoot, 'dependencies'), 'IncludingSubfolders', true));
            testCase.applyFixture(PathFixture(fullfile(repoRoot, 'tests', 'unit', 'utils')));
            folder = testCase.applyFixture(TemporaryFolderFixture).Folder;
            file = fullfile(folder, 'sub_001.mat');
            TestSubproblemRefitTest.writeSubproblem(file, 3);
            testCase.Inputs = load(file);
        end
    end

    methods (Test)
        function testDefault(testCase)
            p = setParams('SILo', struct('isSLAP2', true));
            testCase.verifyEqual(p.debiasVariance, 'fixed');
        end

        function testMissingFieldMeansFixed(testCase)
            d = testCase.Inputs;
            p = d.params;
            if isfield(p, 'debiasVariance'), p = rmfield(p, 'debiasVariance'); end
            out0 = cell(1, 7); out1 = cell(1, 7);
            [out0{:}] = extractSubproblem(d.Y_p, d.F_inv_p, d.sources_p, d.selPix_p, p);
            p.debiasVariance = 'fixed';
            [out1{:}] = extractSubproblem(d.Y_p, d.F_inv_p, d.sources_p, d.selPix_p, p);
            testCase.verifyEqual(out0, out1);
        end

        function testOnlyTheActivityChanges(testCase)
            d = testCase.Inputs;
            p = d.params;
            [Hf, Sf, F0f, LSf, SNRf, ~, Bf] = extractSubproblem(d.Y_p, d.F_inv_p, d.sources_p, d.selPix_p, setfield(p, 'debiasVariance', 'fixed')); %#ok<SFLD>
            [Hm, Sm, F0m, LSm, SNRm, ~, Bm] = extractSubproblem(d.Y_p, d.F_inv_p, d.sources_p, d.selPix_p, setfield(p, 'debiasVariance', 'model')); %#ok<SFLD>
            testCase.verifyEqual(Hf, Hm, 'footprints');
            testCase.verifyEqual(Bf, Bm, 'baseline');
            testCase.verifyEqual(F0f, F0m, 'F0');
            testCase.verifyEqual(LSf, LSm, 'least-squares activity');
            testCase.verifyEqual(Sf ~= 0, Sm ~= 0, 'the same events are refit');
            testCase.verifyLessThan(sum(Sf(:)), sum(Sm(:)), 'smaller activity with the fixed variance');
            testCase.verifyNotEqual(SNRf, SNRm);
        end

        function testSecondChannel(testCase)
            % two-channel data (e.g. spine calcium in channel 2, with its
            % slower kernel params.k2): channel 2's debias refit also holds
            % the variance fixed; its baseline, F0 and least-squares activity
            % do not depend on the option
            d = testCase.Inputs;
            prev = RandStream.setGlobalStream(RandStream('twister', 'Seed', 5));
            restore = onCleanup(@() RandStream.setGlobalStream(prev));
            [r, c] = find(d.selPix_p);
            H = exp(-((r - d.sources_p.R').^2 + (c - d.sources_p.C').^2) / (2*1.4^2)); H = H ./ sum(H, 1);
            spikes = (rand(2, size(d.Y_p, 2)) < 0.01) .* (2 + 4*rand(2, size(d.Y_p, 2))) ./ max(H, [], 1)';
            X2 = convn(spikes, d.params.k2 / max(d.params.k2), 'same');
            Y = cat(3, d.Y_p, d.F_inv_p .* poissrnd((0.6 + H*X2) ./ d.F_inv_p));
            clear restore
            p = d.params;
            outF = cell(1, 7); outM = cell(1, 7);
            [outF{:}] = extractSubproblem(Y, d.F_inv_p, d.sources_p, d.selPix_p, setfield(p, 'debiasVariance', 'fixed')); %#ok<SFLD>
            [outM{:}] = extractSubproblem(Y, d.F_inv_p, d.sources_p, d.selPix_p, setfield(p, 'debiasVariance', 'model')); %#ok<SFLD>
            [Hf, Sf, F0f, LSf, ~, ~, Bf] = outF{:};
            [Hm, Sm, F0m, LSm, ~, ~, Bm] = outM{:};
            testCase.assertEqual(size(Sf, 3), 2);
            testCase.verifyEqual(Hf, Hm, 'footprints');
            testCase.verifyEqual(Bf, Bm, 'baselines, both channels');
            testCase.verifyEqual(F0f, F0m, 'F0, both channels');
            testCase.verifyEqual(LSf, LSm, 'least-squares activity, both channels');
            S2f = Sf(:, :, 2); S2m = Sm(:, :, 2);
            testCase.assertTrue(any(S2m(:) > 0), 'channel 2 has events to refit');
            testCase.verifyEqual(S2f ~= 0, S2m ~= 0, 'channel 2: the same entries are refit');
            testCase.verifyLessThan(sum(S2f(:)), sum(S2m(:)), 'channel 2: smaller activity with the fixed variance');
        end

        function testRejectsUnknownValue(testCase)
            d = testCase.Inputs;
            p = d.params; p.debiasVariance = 'poisson';
            testCase.verifyError(@() extractSubproblem(d.Y_p, d.F_inv_p, d.sources_p, d.selPix_p, p), ...
                'extractTrial:debiasVariance');
            % also in a subproblem without activity, which has nothing to refit
            Yq = d.F_inv_p .* poissrnd(0.8 ./ d.F_inv_p);
            testCase.verifyError(@() extractSubproblem(Yq, d.F_inv_p, d.sources_p, d.selPix_p, p), ...
                'extractTrial:debiasVariance');
        end
    end
end
