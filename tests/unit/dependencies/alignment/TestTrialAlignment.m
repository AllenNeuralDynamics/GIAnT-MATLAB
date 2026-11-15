classdef TestTrialAlignment < matlab.unittest.TestCase
    %TESTTRIALALIGNMENT Unit tests for SILo's across-trial alignment tools:
    %   warpImageAffine, warpTrialImage, trialWarpCoords, trialWarpInverse,
    %   patchFieldMatrix, interpFramesAt, alignTrialsAffine and
    %   alignTrialsPatchwise. The fits are tested on a synthetic labeled
    %   neuron (dendrites with spines, an unscanned band and a hole) warped by
    %   known affine maps and smooth bends, with shot-like noise.
    %   Run with: runtests('tests/unit/dependencies/alignment/TestTrialAlignment.m')

    properties (Constant)
        Sz = [180 260]
        NTrials = 6
    end

    methods (TestClassSetup)
        function addPaths(testCase)
            import matlab.unittest.fixtures.PathFixture
            testsDir = fileparts(mfilename('fullpath'));
            repoRoot = fileparts(fileparts(fileparts(fileparts(testsDir))));
            testCase.applyFixture(PathFixture(fullfile(repoRoot, 'dependencies'), 'IncludingSubfolders', true));
        end
    end

    methods (TestMethodSetup)
        function seedRng(testCase)
            testCase.addTeardown(@rng, rng);
            rng(1, 'twister');
        end
    end

    methods (Test)
        function testRigidWarpMatchesInterp2(testCase)
            % a rigid map resamples exactly as SILo's shift by interp2 did
            IM = rand(50, 70); IM(10:12, 20:30) = nan;
            [cc, rr] = meshgrid(1:70, 1:50); mot = [0.4; -1.3];
            expected = interp2(IM, cc + mot(2), rr + mot(1));
            warp = struct('affine', [1 0 mot(1); 0 1 mot(2)], 'nodes', zeros(0, 2), 'sigma', nan, 'disp', zeros(2, 0), 'patch', []);
            testCase.verifyEqual(warpTrialImage(IM, warp, 1, 'linear'), expected);
            testCase.verifyEqual(isnan(warpTrialImage(IM, warp, 1, 'cubic')), isnan(expected), ...
                'cubic removes no more pixels next to NaN than linear');
        end

        function testFieldReproducesAffineMap(testCase)
            % displacements at the patch centres that follow an affine map give
            % that map everywhere, also far outside the centres
            nodes = [rand(12, 1)*200 + 50, rand(12, 1)*400 + 50];
            A = [1.004 0.002 0.7; -0.003 0.997 -1.2];
            d = (A(:, 1:2)*nodes' + A(:, 3) - nodes')';
            [cc, rr] = meshgrid(linspace(-50, 550, 61), linspace(-50, 350, 41));
            u = patchFieldMatrix(nodes, 30, rr(:), cc(:))*d;
            uA = (A(:, 1:2)*[rr(:) cc(:)]' + A(:, 3) - [rr(:) cc(:)]')';
            testCase.verifyLessThan(max(abs(u - uA), [], 'all'), 1e-9);

            warp = struct('affine', A, 'nodes', nodes, 'sigma', 30, 'disp', d', 'patch', []);
            IM = rand(300, 500);
            a = warpTrialImage(IM, warp, 1, 'linear'); b = warpImageAffine(IM, A, 'linear');
            testCase.verifyEqual(isnan(a), isnan(b));
            testCase.verifyLessThan(max(abs(a(:) - b(:)), [], 'omitnan'), 1e-9);
        end

        function testInverseUndoesBentField(testCase)
            nodes = [rand(12, 1)*200 + 50, rand(12, 1)*400 + 50];
            A = [1.004 0.002 0.7; -0.003 0.997 -1.2];
            d = A(:, 1:2)*nodes' + A(:, 3) - nodes' + 0.4*randn(2, 12);
            warp = struct('affine', A, 'nodes', nodes, 'sigma', 30, 'disp', d, 'patch', []);
            r0 = rand(1000, 1)*300; c0 = rand(1000, 1)*500;
            [rT, cT] = trialWarpCoords(warp, 1, r0, c0);
            [r1, c1] = trialWarpInverse(warp, 1, rT, cT);
            testCase.verifyLessThan(max(abs([r1 - r0; c1 - c0])), 1e-5);
        end

        function testInterpFramesAtMatchesInterpFrames(testCase)
            % on a separable grid (a shift), the point-wise interpolation is
            % the separable one
            M1 = rand(60, 80, 2, 'single')*100; M1(1:5, :, :) = nan; M1(30, 40, :) = nan;
            F = rand(60, 80, 'single') + 0.5; F(isnan(M1(:, :, 1))) = 0;
            for sh = [0 0; 0.3 -0.7; -2.25 1.5]'
                viewC = (1:80) - 4 + sh(2); viewR = (1:60)' - 4 + sh(1);
                [a, va] = interpFrames(M1, viewC, viewR, F);
                [cg, rg] = meshgrid(viewC, viewR);
                [b, vb] = interpFramesAt(M1, rg, cg, F);
                b = reshape(b, 60, 80, 2); vb = reshape(vb, 60, 80);
                testCase.verifyEqual(isnan(b), isnan(a));
                testCase.verifyLessThan(max(abs(a(:) - b(:)), [], 'omitnan'), 1e-3);
                testCase.verifyLessThan(max(abs(va(:) - vb(:)), [], 'omitnan'), 1e-4);
            end
        end

        function testAffineRecoversKnownMaps(testCase)
            [truth, scan] = TestTrialAlignment.neuron();
            [M, ~, U, A0] = TestTrialAlignment.trials(truth, scan, false);
            A = alignTrialsAffine(M, A0);
            [err, err0] = TestTrialAlignment.fieldError(@(t, r, c) affineCoords(A(:, :, t), r, c), U, A0, truth, scan);
            testCase.verifyLessThan(err, 0.05, 'relative maps between trials recovered to 0.05 px');
            testCase.verifyLessThan(err, err0/5, 'much better than the rigid start');
        end

        function testPatchwiseRecoversBends(testCase)
            % smooth bends (over ~40 px) on top of affine maps
            [truth, scan] = TestTrialAlignment.neuron();
            [M, Mh, U, A0] = TestTrialAlignment.trials(truth, scan, true);
            [warp, info] = alignTrialsPatchwise(M, A0, Mh, 1500);
            [err, err0] = TestTrialAlignment.fieldError(@(t, r, c) trialWarpCoords(warp, t, r, c), U, A0, truth, scan);
            testCase.verifyLessThan(err, 0.08, 'relative fields between trials recovered to 0.08 px');
            testCase.verifyLessThan(err, err0/3, 'much better than the rigid start');
            testCase.verifyGreaterThan(size(warp.nodes, 1), 10);
            testCase.verifyGreaterThanOrEqual(min(info.nPx), 750, 'no patch under half the target size');
            Aonly = alignTrialsAffine(M, A0);
            errA = TestTrialAlignment.fieldError(@(t, r, c) affineCoords(Aonly(:, :, t), r, c), U, A0, truth, scan);
            testCase.verifyLessThan(err, 0.7*errA, 'the bends are captured beyond an affine map');
        end

        function testPatchwiseAddsNothingWithoutBends(testCase)
            % affine maps only: the patchwise field is as good as the affine fit
            [truth, scan] = TestTrialAlignment.neuron();
            [M, Mh, U, A0] = TestTrialAlignment.trials(truth, scan, false);
            warp = alignTrialsPatchwise(M, A0, Mh, 1500);
            err = TestTrialAlignment.fieldError(@(t, r, c) trialWarpCoords(warp, t, r, c), U, A0, truth, scan);
            testCase.verifyLessThan(err, 0.05);
        end
    end

    methods (Static)
        function [truth, scan] = neuron()
            % mean image of a labeled neuron: background, dendrites as smooth
            % curves and spines as blobs next to them; an unscanned band and a
            % hole are NaN
            sz = TestTrialAlignment.Sz;
            [cc, rr] = meshgrid(1:sz(2), 1:sz(1));
            truth = 5*ones(sz);
            for k = 1:6
                t = linspace(0, 1, 400)';
                p0 = rand(1, 2).*sz; ang = rand*pi; len = 150 + 100*rand;
                path = p0 + t*len*[sin(ang) cos(ang)] + 8*sin(2*pi*(t*(1 + rand) + rand))*[cos(ang) -sin(ang)];
                for q = 1:20:size(path, 1)
                    truth = truth + 30*exp(-((rr - path(q, 1)).^2 + (cc - path(q, 2)).^2)/(2*1.5^2));
                    sp = path(q, :) + (3 + 2*rand)*[cos(ang) -sin(ang)]*sign(randn);
                    truth = truth + (40 + 80*rand)*exp(-((rr - sp(1)).^2 + (cc - sp(2)).^2)/(2*1.1^2));
                end
            end
            scan = true(sz); scan(:, 1:12) = false; scan(80:100, 120:150) = false;
            truth(~scan) = nan;
        end

        function [M, Mh, U, A0] = trials(truth, scan, bends)
            % trial mean images (and two noisy halves of each) of the truth seen
            % through known maps: common-frame pixel x is at x + U(x) in the trial
            sz = size(truth); nT = TestTrialAlignment.NTrials;
            [cc, rr] = meshgrid(1:sz(2), 1:sz(1));
            filled = truth; filled(isnan(filled)) = 5;
            M = nan([sz nT]); Mh = nan([sz nT 2]); U = nan([sz 2 nT]); A0 = zeros(2, 3, nT);
            for t = 1:nT
                L = eye(2) + 0.003*randn(2); tr = 0.6*randn(2, 1);
                uR = (L(1, 1) - 1)*rr + L(1, 2)*cc + tr(1);
                uC = L(2, 1)*rr + (L(2, 2) - 1)*cc + tr(2);
                if bends
                    bR = imgaussfilt(randn(sz), 40, 'FilterSize', 401); bC = imgaussfilt(randn(sz), 40, 'FilterSize', 401);
                    uR = uR + 0.25*bR/std(bR(:)); uC = uC + 0.25*bC/std(bC(:));
                end
                U(:, :, :, t) = cat(3, uR, uC);
                FR = griddedInterpolant(uR, 'linear', 'nearest'); FC = griddedInterpolant(uC, 'linear', 'nearest');
                y = [rr(:) cc(:)]';
                for it = 1:8 % the common-frame point each trial pixel shows
                    y = [rr(:)'; cc(:)'] - [FR(y(1, :), y(2, :)); FC(y(1, :), y(2, :))];
                end
                img = reshape(interp2(filled, y(2, :), y(1, :), 'cubic', 5), sz);
                for h = 1:2
                    x = img + 0.2*sqrt(max(img, 0)).*randn(sz); x(~scan) = nan; % scan mask fixed in the trial's coordinates
                    Mh(:, :, t, h) = x;
                end
                M(:, :, t) = mean(Mh(:, :, t, :), 4);
                c0 = round(sz/2);
                A0(:, :, t) = [1 0 round(uR(c0(1), c0(2))); 0 1 round(uC(c0(1), c0(2)))];
            end
        end

        function [err, err0] = fieldError(coords, U, A0, truth, scan)
            % RMS error of the fitted maps relative to the trials' mean map
            % (the common frame is the trials' average geometry), weighted by
            % image structure (squared gradient of the sqrt image), which is
            % where the data constrain the map; and the same for the rigid start
            sz = size(truth); nT = size(U, 4);
            [cc, rr] = meshgrid(1:sz(2), 1:sz(1));
            sel = imerode(scan, ones(5));
            est = nan(nnz(sel), 2, nT); tru = est; rig = est;
            for t = 1:nT
                [rT, cT] = coords(t, rr(sel), cc(sel));
                est(:, :, t) = [rT - rr(sel), cT - cc(sel)];
                uR = U(:, :, 1, t); uC = U(:, :, 2, t);
                tru(:, :, t) = [uR(sel), uC(sel)];
                rig(:, :, t) = repmat(A0(:, 3, t)', nnz(sel), 1);
            end
            rel = @(x) x - mean(x, 3);
            [gc, gr] = gradient(sqrt(max(truth, 0))); w = gr(sel).^2 + gc(sel).^2; w(isnan(w)) = 0;
            wrms = @(e) sqrt(sum(mean(sum(e.^2, 2), 3).*w)/sum(w));
            err = wrms(rel(est) - rel(tru));
            err0 = wrms(rel(rig) - rel(tru));
        end
    end
end

function [rT, cT] = affineCoords(A, r, c)
rT = A(1, 1)*r + A(1, 2)*c + A(1, 3);
cT = A(2, 1)*r + A(2, 2)*c + A(2, 3);
end
