classdef TestSubproblemSaving < matlab.unittest.TestCase
    %TESTSUBPROBLEMSAVING extractTrial's developer option
    %params.saveSubproblemsDir: with it set, the inputs of every
    %subproblem are saved (sub_###.mat, one folder per trial), and
    %extractSubproblem on a saved file reproduces that subproblem's part of
    %the trial's result exactly. Saving does not change the result, and it
    %is off by default. Uses the real extractTrial on a small synthetic
    %trial (two subproblems), on a one-worker process pool.
    %   Run with: runtests('tests/unit/dependencies/extraction/TestSubproblemSaving.m')

    properties
        Folder
    end

    methods (TestClassSetup)
        function addPathsAndPool(testCase)
            import matlab.unittest.fixtures.PathFixture
            testsDir = fileparts(mfilename('fullpath'));
            repoRoot = fileparts(fileparts(fileparts(fileparts(testsDir))));
            testCase.applyFixture(PathFixture(fullfile(repoRoot, 'dependencies'), 'IncludingSubfolders', true));
            if isempty(gcp('nocreate'))
                extractionPool(struct('nWorkers', 1));
                testCase.addTeardown(@() delete(gcp('nocreate')));
            end
        end
    end

    methods (TestMethodSetup)
        function makeFolder(testCase)
            import matlab.unittest.fixtures.TemporaryFolderFixture
            testCase.Folder = testCase.applyFixture(TemporaryFolderFixture).Folder;
        end
    end

    methods (Test)
        function testOffByDefault(testCase)
            p = setParams('SILo', struct('isSLAP2', true));
            testCase.verifyEqual(p.saveSubproblemsDir, '');
            testCase.verifyEqual(subproblemSaveDir(p), '');
        end

        function testSaveDirValues(testCase)
            % off: missing, empty, blank, false or 0 (e.g. from JSON); text: on
            for v = {'', "", '  ', [], false, 0}
                testCase.verifyEqual(subproblemSaveDir(struct('saveSubproblemsDir', v)), '');
            end
            testCase.verifyEqual(subproblemSaveDir(struct()), '');
            testCase.verifyEqual(subproblemSaveDir(struct('saveSubproblemsDir', " D:\dumps ")), 'D:\dumps');
            testCase.verifyError(@() subproblemSaveDir(struct('saveSubproblemsDir', true)), 'extractTrial:saveSubproblemsDir');
            testCase.verifyError(@() subproblemSaveDir(struct('saveSubproblemsDir', 3)), 'extractTrial:saveSubproblemsDir');
        end

        function testSavedInputsReproduceTheTrial(testCase)
            [Y, Finv, sources, selPix, params] = TestSubproblemSaving.syntheticTrial();
            params.saveSubproblemsDir = testCase.Folder;
            params.subproblemSource = struct('path', 2, 'trial', 3, 'file', fullfile('raw', 'acq_DMD2-TRIAL000004.dat'));
            [H, S, LS, F0, SNR] = fetchOutputs(extractTrial(Y, Finv, sources, selPix, params));

            d = fullfile(testCase.Folder, 'path2_trial003_acq_DMD2-TRIAL000004');
            testCase.assertTrue(isfolder(d), 'one folder per trial, named by path, trial and file');
            t = load(fullfile(d, 'trial.mat'));
            testCase.verifyEqual(t.sources, sources);
            testCase.verifyEqual(t.selPix, selPix);
            testCase.verifyEqual(t.nProblems, 2);
            files = dir(fullfile(d, 'sub_*.mat'));
            testCase.assertEqual({files.name}, {'sub_001.mat', 'sub_002.mat'}, 'one file per subproblem');

            selIdxs = zeros(size(selPix)); selIdxs(selPix) = 1:nnz(selPix);
            seen = [];
            for k = 1:numel(files)
                s = load(fullfile(d, files(k).name));
                testCase.verifyEqual(s.problemIx, k);
                testCase.verifyEqual(s.nProblems, 2);
                testCase.verifyEqual(s.params.photonScale, params.photonScale);
                testCase.verifyEqual([s.srcR s.srcC], [sources.R(s.srcIdx) sources.C(s.srcIdx)]);
                testCase.verifyEqual([s.sources_p.R s.sources_p.C], ...
                    [s.srcR - s.rowSupport(1) + 1, s.srcC - s.colSupport(1) + 1], 'cropped source positions');
                full = false(size(selPix));
                full(s.rowSupport(1):s.rowSupport(2), s.colSupport(1):s.colSupport(2)) = s.selPix_p;
                rows = selIdxs(full);
                testCase.verifyEqual(s.Y_p, Y(rows, :) / params.photonScale, 'photon-normalized data, as fitted');
                testCase.verifyEqual(s.F_inv_p, Finv(rows, :));
                % solved on a pool worker, as extractTrial solves it
                [Hi, Si, F0i, LSi, SNRi] = fetchOutputs(parfeval(@extractSubproblem, 5, ...
                    s.Y_p, s.F_inv_p, s.sources_p, s.selPix_p, s.params));
                testCase.verifyEqual(H(rows, s.srcIdx), single(Hi), 'footprints');
                testCase.verifyEqual(S(s.srcIdx, :, :), Si .* params.photonScale, 'events');
                testCase.verifyEqual(LS(s.srcIdx, :, :), LSi .* params.photonScale, 'least-squares activity');
                testCase.verifyEqual(F0(s.srcIdx, :, :), F0i .* params.photonScale, 'F0');
                testCase.verifyEqual(SNR(s.srcIdx), SNRi, 'SNR');
                seen = [seen; s.srcIdx(:)]; %#ok<AGROW>
            end
            testCase.verifyEqual(sort(seen), (1:numel(sources.R))', 'every source in exactly one subproblem');
        end

        function testSavingLeavesTheResultUnchanged(testCase)
            [Y, Finv, sources, selPix, params] = TestSubproblemSaving.syntheticTrial();
            out0 = cell(1, 5); out1 = cell(1, 5);
            [out0{:}] = fetchOutputs(extractTrial(Y, Finv, sources, selPix, params));
            params.saveSubproblemsDir = testCase.Folder;
            [out1{:}] = fetchOutputs(extractTrial(Y, Finv, sources, selPix, params));
            testCase.verifyEqual(out1, out0);
        end

        function testSavesGroundTruth(testCase)
            % with ground truth given (extractTrial's 7th argument), each
            % subproblem's slice of it is saved too, in the solver's units
            [Y, Finv, sources, selPix, params, GT] = TestSubproblemSaving.syntheticTrial();
            params.saveSubproblemsDir = testCase.Folder;
            set(0, 'DefaultFigureVisible', 'off'); % the ground-truth path plots its progress
            testCase.addTeardown(@() set(0, 'DefaultFigureVisible', 'on'));
            fetchOutputs(extractTrial(Y, Finv, sources, selPix, params, [], GT));
            d = dir(fullfile(testCase.Folder, '*', 'sub_*.mat'));
            testCase.assertNumElements(d, 2);
            selIdxs = zeros(size(selPix)); selIdxs(selPix) = 1:nnz(selPix);
            ps = params.photonScale;
            for k = 1:numel(d)
                s = load(fullfile(d(k).folder, d(k).name));
                full = false(size(selPix));
                full(s.rowSupport(1):s.rowSupport(2), s.colSupport(1):s.colSupport(2)) = s.selPix_p;
                rows = selIdxs(full);
                testCase.verifyEqual(s.GTp.S, GT.S(s.srcIdx, :) / ps, 'events, photon-normalized');
                testCase.verifyEqual(s.GTp.X, GT.X(s.srcIdx, :) / ps);
                testCase.verifyEqual(s.GTp.B, GT.B(rows, :) / ps, 'baseline at the subproblem''s pixels');
                testCase.verifyEqual(s.GTp.H, GT.H(rows, s.srcIdx), 'footprints at its pixels and sources');
                testCase.verifyEqual(s.GTp.Hs, GT.Hs(rows, s.srcIdx));
            end
        end

        function testUnnamedTrialsDoNotCollide(testCase)
            % without params.subproblemSource (extractTrial called directly),
            % each call gets a folder of its own
            [Y, Finv, sources, selPix, params] = TestSubproblemSaving.syntheticTrial();
            params.saveSubproblemsDir = testCase.Folder;
            fetchOutputs(extractTrial(Y, Finv, sources, selPix, params));
            fetchOutputs(extractTrial(Y, Finv, sources, selPix, params));
            d = dir(testCase.Folder); d = d([d.isdir] & ~startsWith({d.name}, '.'));
            testCase.verifyNumElements(d, 2);
            for k = 1:numel(d)
                testCase.verifyNumElements(dir(fullfile(testCase.Folder, d(k).name, 'sub_*.mat')), 2);
            end
        end
    end

    methods (Static)
        function [Y, Finv, sources, selPix, params, GT] = syntheticTrial()
            % three sources in two separate groups, low photon counts (about
            % one photon per pixel and sample, as on SLAP2) and sparse events
            prev = RandStream.setGlobalStream(RandStream('twister', 'Seed', 7)); % poissrnd draws from the global stream
            restore = onCleanup(@() RandStream.setGlobalStream(prev));
            sz = [32 32];
            sources = struct('R', [9; 11; 23], 'C', [9; 12; 21], 'V', ones(3, 1));
            params = setParams('SILo', struct('isSLAP2', true, 'nWorkers', 1));
            params.minBaseline = 0.01; % what SILo sets before extraction
            params.numChannels = 1;
            params.selRadius = ceil(2*params.dXY);
            params.tau_full = params.tau_s*params.analyzeHz;
            params = setParamsExtractTrial(params);
            selPix = false(sz); selPix(sub2ind(sz, sources.R, sources.C)) = true;
            selPix = imdilate(selPix, strel('disk', params.selRadius));
            [r, c] = find(selPix);
            nT = 600;
            H = exp(-((r - sources.R').^2 + (c - sources.C').^2) / (2*1.4^2));
            H = H ./ sum(H, 1);
            spikes = (rand(3, nT) < 0.02) .* (20 + 30*rand(3, nT));
            X = convn(spikes, params.k / max(params.k), 'same');
            B = 0.8 + 0.2*sin((1:nT)/90) + 0.1*rand(numel(r), 1);
            Finv = 0.77*ones(numel(r), nT);
            T = B + H*X;
            Y = params.photonScale * Finv .* poissrnd(T ./ Finv); % counts; extractTrial divides by photonScale
            % the ground truth, in counts as extractTrial takes it
            ps = params.photonScale;
            GT = struct('S', ps * spikes / max(params.k), 'X', ps * X, 'B', ps * B, 'Y', Y, 'H', H, 'Hs', H);
        end
    end
end
