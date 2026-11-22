classdef TestExtractTrialUnmeasured < matlab.unittest.TestCase
    %TESTEXTRACTTRIALUNMEASURED Source extraction (extractTrial) when a
    %source's pixels have no samples in part or all of a trial, as at the
    %edge of the field of view after motion or an online offset change. A
    %small simulated recording: three sources on a dim background with a
    %bright dendrite, Poisson noise; the edge source's pixels (columns 20 on)
    %are missing in a block of frames or in the whole trial, and are filled
    %in as processAllTrials_Async's trial loader does (a moving mean of the
    %pixel's samples, or 0, with 1000 times the mean variance multiplier).
    %The fit must not estimate the edge source's activity where its
    %footprint is mostly unmeasured, and must leave the other sources as
    %they were. Runs extractTrial's subproblems on a thread pool.
    %   Run with: runtests('tests/unit/dependencies/extraction/TestExtractTrialUnmeasured.m')

    properties (Constant)
        Size = [31 27]
        NumFrames = 3742
        Sources = [16 21; 16 12; 9 15] % [row col]; the first is at the edge
        EdgeCol = 20                     % columns >= EdgeCol can be missing
        Block = 2001:3000                % frames missing in 'block' mode
    end

    properties
        Params
    end

    methods (TestClassSetup)
        function addPaths(testCase)
            import matlab.unittest.fixtures.PathFixture
            testsDir = fileparts(mfilename('fullpath'));
            repoRoot = fileparts(fileparts(fileparts(fileparts(testsDir))));
            testCase.applyFixture(PathFixture(fullfile(repoRoot, 'dependencies'), 'IncludingSubfolders', true));
        end

        function requireToolboxes(testCase)
            testCase.assumeNotEmpty(ver('images'), 'Requires Image Processing Toolbox.');
            testCase.assumeNotEmpty(ver('optim'), 'Requires Optimization Toolbox.');
            testCase.assumeNotEmpty(ver('parallel'), 'Requires Parallel Computing Toolbox.');
        end

        function startPool(testCase)
            pool = gcp('nocreate');
            if isempty(pool)
                parpool('Threads', 2);
                testCase.addTeardown(@() delete(gcp('nocreate')));
            end
        end

        function makeParams(testCase)
            params = setParams('SILo', struct('isSLAP2', true, 'nWorkers', 2));
            params.minBaseline = 0.01; % as SILo
            params.selRadius = ceil(2*params.dXY);
            params.tau_full = params.tau_s*params.analyzeHz;
            testCase.Params = setParamsExtractTrial(params);
        end
    end

    methods (Test)
        function testBlockOfUnmeasuredFrames(testCase)
            [Y, Finv, measured, sources, selPix] = testCase.recording('block');
            [~, S, LS, F0, SNR] = testCase.extract(Y, Finv, sources, selPix, measured);
            missing = false(1, testCase.NumFrames); missing(testCase.Block) = true;
            testCase.verifyEqual(isnan(S(1,:)), missing, 'edge source: activity not estimated exactly where it is unmeasured');
            testCase.verifyEqual(isnan(LS(1,:)), missing);
            testCase.verifyEqual(isnan(F0(1,:)), missing);
            testCase.verifyTrue(isfinite(SNR(1)), 'edge source: SNR from the frames it is measured in');
            testCase.verifyTrue(all(isfinite(S(2:3,:)), 'all') && all(isfinite(LS(2:3,:)), 'all') && all(isfinite(F0(2:3,:)), 'all'), ...
                'sources inside the field are estimated in every frame');
        end

        function testSourceUnmeasuredInWholeTrial(testCase)
            [Y, Finv, measured, sources, selPix] = testCase.recording('trial');
            [~, S, LS, F0, SNR] = testCase.extract(Y, Finv, sources, selPix, measured);
            testCase.verifyTrue(all(isnan(S(1,:))) && all(isnan(LS(1,:))) && all(isnan(F0(1,:))) && isnan(SNR(1)), ...
                'edge source: not estimated in a trial that never measures it');
            testCase.verifyTrue(all(isfinite(S(2:3,:)), 'all') && all(isfinite(SNR(2:3))));
        end

        function testFullyMeasuredUnchanged(testCase)
            % with every sample measured, passing the mask changes nothing
            [Y, Finv, measured, sources, selPix] = testCase.recording('none');
            testCase.assertTrue(all(measured(:)));
            [H1, S1, LS1, F01, SNR1] = testCase.extract(Y, Finv, sources, selPix, measured);
            [H0, S0, LS0, F00, SNR0] = testCase.extract(Y, Finv, sources, selPix, []);
            testCase.verifyEqual({H1, S1, LS1, F01, SNR1}, {H0, S0, LS0, F00, SNR0});
        end
    end

    methods
        function [Y, Finv, measured, sources, selPix] = recording(testCase, mode)
            % Y: counts at the selPix pixels (pixels x frames), filled in where
            % missing, as the trial loader hands them to extractTrial
            rng(1, 'twister');
            p = testCase.Params; sz = testCase.Size; T = testCase.NumFrames; src = testCase.Sources;
            [cc, rr] = meshgrid(1:sz(2), 1:sz(1));
            imaged = cc <= 23; % the canvas' imaged region
            selPix = false(sz);
            selPix(sub2ind(sz, src(:,1), src(:,2))) = true;
            selPix = imdilate(selPix, strel('disk', p.selRadius)) & imaged;
            background = 0.15 + 0.5*exp(-(rr - 16).^2/(2*2^2)); % photons per sample
            H = zeros(prod(sz), size(src, 1));
            for s = 1:size(src, 1)
                h = exp(-((rr - src(s,1)).^2 + (cc - src(s,2)).^2)/(2*1.6^2));
                H(:, s) = h(:)/sum(h(:));
            end
            events = (rand(size(src, 1), T) < 2/p.analyzeHz) .* (50 + 100*rand(size(src, 1), T));
            X = convn(events, p.k, 'same');
            Y = poissrnd(background(:) + H*X) * p.photonScale;
            Y = Y(selPix(:), :);
            Finv = ones(size(Y));
            missing = false(size(Y));
            atEdge = cc(selPix) >= testCase.EdgeCol;
            switch mode
                case 'block', missing(atEdge, testCase.Block) = true;
                case 'trial', missing(atEdge, :) = true;
            end
            measured = ~missing;
            % fill in as processAllTrials_Async's loadTrial
            Y(missing) = nan;
            fill = Y; fill(all(missing, 2), :) = 0;
            while any(isnan(fill), 'all')
                fill = smoothdata(fill, 2, 'movmean', p.baselineWindow_samps, 'omitmissing');
            end
            Y(missing) = fill(missing);
            Finv(missing) = 1000*mean(Finv, 'all');
            sources = struct('R', src(:,1), 'C', src(:,2));
        end

        function [H, S, LS, F0, SNR] = extract(testCase, Y, Finv, sources, selPix, measured)
            if isempty(measured)
                fut = extractTrial(Y, Finv, sources, selPix, testCase.Params);
            else
                fut = extractTrial(Y, Finv, sources, selPix, testCase.Params, measured);
            end
            [H, S, LS, F0, SNR] = fetchOutputs(fut);
        end
    end
end
