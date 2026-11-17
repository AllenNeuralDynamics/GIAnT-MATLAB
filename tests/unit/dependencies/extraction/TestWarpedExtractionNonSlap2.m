classdef TestWarpedExtractionNonSlap2 < matlab.unittest.TestCase
    %TESTWARPEDEXTRACTIONNONSLAP2 processAllTrials_Async's loader for
    %non-SLAP2 movies (an h5 of rows x cols x channels x frames) with an
    %affine or bent patchwise map from the common frame to the trial. Each
    %selected common-frame pixel should be the bilinear interpolation of each
    %baseline-subtracted frame at its mapped position, with the channels in
    %analysis order. A sample is missing only where a mapped neighbour is
    %NaN, not where the unwarped frame has a NaN at the pixel's own index.
    %The samples are those handed to source extraction (extractTrial and
    %extractionPool are test doubles in ./fakes).
    %   Run with: runtests('tests/unit/dependencies/extraction/TestWarpedExtractionNonSlap2.m')

    properties (Constant)
        Sz = [40 40]
        NumFrames = 120 % processAllTrials_Async extracts only trials of at least 100 frames
        Hole = {20:21, 18:19} % rows, cols never imaged
    end

    properties
        DataDir
        Name = 'trial1_REGISTERED.h5'
    end

    methods (TestClassSetup)
        function addPaths(testCase)
            import matlab.unittest.fixtures.PathFixture
            testsDir = fileparts(mfilename('fullpath'));
            repoRoot = fileparts(fileparts(fileparts(fileparts(testsDir))));
            testCase.applyFixture(PathFixture(fullfile(repoRoot, 'dependencies'), 'IncludingSubfolders', true));
            % added last so that the doubles shadow the real extractor and pool
            testCase.applyFixture(PathFixture(fullfile(testsDir, 'fakes')));
        end
    end

    methods (TestMethodSetup)
        function makeFolder(testCase)
            import matlab.unittest.fixtures.TemporaryFolderFixture
            testCase.DataDir = testCase.applyFixture(TemporaryFolderFixture).Folder;
            testCase.addTeardown(@() rmappdata0('FakeExtractTrialFile'));
        end
    end

    methods (Test)
        function testSamplesAtMappedPositions(testCase)
            [mov, tex, gain] = TestWarpedExtractionNonSlap2.movie();
            h5create(fullfile(testCase.DataDir, testCase.Name), '/data', size(mov));
            h5write(fullfile(testCase.DataDir, testCase.Name), '/data', mov);
            base = prctile(reshape(median(mov, 4, 'omitmissing'), [], 2), 5, 1); % the loader's baseline per channel
            selPix = false(TestWarpedExtractionNonSlap2.Sz); selPix(6:35, 6:35) = true;
            [r, c] = find(selPix);
            inHole = ismember(r, TestWarpedExtractionNonSlap2.Hole{1}) & ismember(c, TestWarpedExtractionNonSlap2.Hole{2});
            order = [2 1]; % activityChannel 2: analysis order puts channel 2 first
            maps = {TestWarpedExtractionNonSlap2.affineWarp(), TestWarpedExtractionNonSlap2.bentWarp()};
            names = {'affine', 'patchwise'};
            for w = 1:2
                [Y, Finv] = testCase.extract(selPix, maps{w});
                if w == 1 % affine: the mapped positions straight from the map
                    A = maps{w}.affine; P = A(:, 1:2)*[r c]' + A(:, 3); rT = P(1, :)'; cT = P(2, :)';
                else % patchwise: the field itself is tested in TestTrialAlignment
                    [rT, cT] = trialWarpCoords(maps{w}, 1, r, c);
                end
                missing = isnan(interp2(tex(:, :, 1), cT, rT, 'linear'));
                testCase.verifyEqual(nnz(missing), nnz(any(isnan(interp2(tex(:, :, 2), cT, rT, 'linear')), 2)), names{w});
                testCase.assertGreaterThan(nnz(missing), 0, [names{w} ': the hole is sampled']);
                filled = Finv > 100*median(Finv(:)); % samples the loader found NaN (it fills them and inflates their variance)
                testCase.verifyEqual(all(filled, 2), missing, [names{w} ': missing exactly where a mapped neighbour is NaN']);
                testCase.verifyEqual(any(filled, 2), missing, [names{w} ': in every frame']);
                testCase.assertTrue(any(inHole & ~missing), [names{w} ': some pixels at the hole''s indices map away from it']);
                testCase.verifyFalse(any(filled(inHole & ~missing, :), 'all'), [names{w} ': those pixels are sampled']);
                testCase.verifyEqual(Finv(~missing, :), ones(nnz(~missing), TestWarpedExtractionNonSlap2.NumFrames), 'AbsTol', 1e-6, ... % the loader keeps single precision
                    [names{w} ': variance multiplier 1 for unit freshness']);
                for k = 1:2
                    expected = interp2(tex(:, :, order(k)), cT, rT, 'linear') .* gain - base(order(k));
                    testCase.verifyEqual(Y(~missing, :, k), expected(~missing, :), 'AbsTol', 1e-5, ...
                        sprintf('%s: channel %d sampled at the mapped positions', names{w}, order(k)));
                end
            end
        end
    end

    methods (Static)
        function [mov, tex, gain] = movie()
            % two smooth static textures scaled by a gain that varies over
            % frames (the loader drops samples far below their temporal SD,
            % so a constant movie would lose its lowest pixels), NaN in the hole
            sz = TestWarpedExtractionNonSlap2.Sz; nFr = TestWarpedExtractionNonSlap2.NumFrames;
            [C, R] = meshgrid(1:sz(2), 1:sz(1));
            tex = cat(3, 100 + 10*sin(R/3).*cos(C/4), 100 + 8*cos(R/5 + C/7));
            tex(TestWarpedExtractionNonSlap2.Hole{:}, :) = nan;
            gain = 1 + 0.2*sin((1:nFr)/5);
            mov = tex .* reshape(gain, 1, 1, 1, nFr);
        end

        function warp = affineWarp()
            % a stretch, shear and shift of about 2 px, so that the hole's own
            % pixels map away from it
            A = [1.03 0.02; -0.015 0.97]; x0 = [20.5; 20.5];
            warp = struct('affine', [A, x0 - A*x0 + [2.4; -1.7]], 'nodes', zeros(0, 2), 'sigma', nan, 'disp', zeros(2, 0), 'patch', []);
        end

        function warp = bentWarp()
            % the affine map plus bends of up to 0.6 px at patch centres
            a = TestWarpedExtractionNonSlap2.affineWarp();
            [nc, nr] = meshgrid([14 20 27], [13 20 28]); nodes = [nr(:) nc(:)];
            bend = 0.6*[1 -0.5 0 0.3 -1 0.4 0 0.8 -0.6; -0.4 0.7 -1 0 0.5 0 0.9 -0.3 0.2];
            d = a.affine(:, 1:2)*nodes' + a.affine(:, 3) - nodes' + bend;
            warp = struct('affine', a.affine, 'nodes', nodes, 'sigma', 3, 'disp', d, 'patch', []);
        end
    end

    methods
        function [Y, Finv] = extract(testCase, selPix, warp)
            [~, tex] = TestWarpedExtractionNonSlap2.movie();
            [r, c] = find(selPix);
            sources = struct('R', r, 'C', c, 'V', ones(size(r)));
            params = setParams('SILo', struct('isSLAP2', false, 'nWorkers', 1, 'activityChannel', 2));
            % the parameters SILo derives before extraction
            params.minBaseline = 0.01;
            params.numChannels = 2;
            params.alignHz = 20;
            params.selRadius = ceil(2 * params.dXY);
            params.tau_full = params.tau_s * params.analyzeHz;
            params = setParamsExtractTrial(params);
            out = [tempname '.mat'];
            testCase.addTeardown(@() deleteIfExists(out));
            setappdata(0, 'FakeExtractTrialFile', out);
            nFr = TestWarpedExtractionNonSlap2.NumFrames;
            alignData = struct('motionC', zeros(1, nFr), 'motionR', zeros(1, nFr));
            processAllTrials_Async(testCase.DataDir, {testCase.Name}, 1, nFr, selPix, sources, ...
                {false}, {alignData}, median(tex, 4), warp, {}, 1, params);
            S = load(out);
            Y = S.Y_obs;
            Finv = S.Finv;
        end
    end
end

function rmappdata0(name)
if isappdata(0, name), rmappdata(0, name); end
end

function deleteIfExists(fn)
if isfile(fn), delete(fn); end
end
