classdef TestFindSourcesWMF < matlab.unittest.TestCase
    %TESTFINDSOURCESWMF Unit tests for findSourcesWMF (params.peakMethod 'wmf').
    %   findSourcesWMF seeds sources from the trials' aligned activity images:
    %   a whitened matched filter of their combined robust z image, with
    %   nearby maxima merged unless they are resolved in both halves of the
    %   trials, seeds added where the residual of their fit still peaks in
    %   both halves, and seeds kept only if active in at least 2 trials. These
    %   tests use synthetic activity images: noise correlated along rows (as
    %   along scan lines), Gaussian sources whose amplitude varies by trial,
    %   and NaN where pixels were not imaged.
    %   Run with: runtests('tests/unit/dependencies/TestFindSourcesWMF.m')

    properties (Constant)
        Sz = [96 96]
        NTrials = 10
        Amp = 4 % mean source peak, in units of the per-trial noise SD (weaker sources are often too faint in single trials to pass the 2-trial activity check)
    end

    properties (TestParameter)
        nFew = struct('oneTrial', 1, 'twoTrials', 2, 'threeTrials', 3)
    end

    methods (TestClassSetup)
        function addPaths(testCase)
            import matlab.unittest.fixtures.PathFixture
            testsDir = fileparts(mfilename('fullpath'));
            repoRoot = fileparts(fileparts(fileparts(testsDir)));
            testCase.applyFixture(PathFixture(fullfile(repoRoot, 'dependencies'), 'IncludingSubfolders', true));
        end
    end

    methods (TestMethodSetup)
        function seedRng(testCase)
            testCase.addTeardown(@rng, rng);
            rng(0, 'twister');
        end
    end

    methods (Test)
        function testFindsIsolatedSources(testCase)
            src = TestFindSourcesWMF.grid();
            act = TestFindSourcesWMF.simulate(src, ones(size(src, 1), 1), TestFindSourcesWMF.NTrials);
            [seeds, info] = findSourcesWMF(act, TestFindSourcesWMF.params());
            d = TestFindSourcesWMF.nearest(src, seeds);
            testCase.verifyLessThan(max(d), 1, 'every source has a seed within 1 px');
            testCase.verifyEqual(size(seeds, 1), size(src, 1), 'one seed per source and no extras');
            testCase.verifyEqual(size(info.F), TestFindSourcesWMF.Sz);
            testCase.verifyTrue(all(diff(info.value) <= 0), 'seeds are ordered strongest first');
        end

        function testElongatedSourceGivesOneSeed(testCase)
            % a source 2.5x longer than wide, which a Gaussian fit can split
            % into seeds at its ends, shows no dip in the filtered image, so
            % before deblending it is one seed at its centre
            [act, c] = TestFindSourcesWMF.withElongatedSource();
            p = TestFindSourcesWMF.params(); p.wmfDeblendth = Inf;
            seeds = findSourcesWMF(act, p);
            near = hypot(seeds(:, 1) - c(1), seeds(:, 2) - c(2)) < 6;
            testCase.verifyEqual(nnz(near), 1, 'one seed on the elongated source');
            testCase.verifyLessThan(hypot(seeds(near, 1) - c(1), seeds(near, 2) - c(2)), 1, 'at its centre');
        end

        function testDeblendingKeepsSeedsOnElongatedSource(testCase)
            % deblending fits the source with copies of the shape of the
            % others, and its residual reproduces in both halves, so a source
            % this much longer than the others can get a second seed; seeds
            % stay on it, one of them near its centre
            [act, c] = TestFindSourcesWMF.withElongatedSource();
            [seeds, info] = findSourcesWMF(act, TestFindSourcesWMF.params());
            near = hypot(seeds(:, 1) - c(1), seeds(:, 2) - c(2)) < 6;
            testCase.verifyLessThanOrEqual(nnz(near), 2, 'at most two seeds on the elongated source');
            testCase.verifyLessThan(abs(seeds(near, 1) - c(1)), 1, 'on its axis');
            testCase.verifyLessThan(abs(seeds(near, 2) - c(2)), 3, 'within one length scale of its centre');
            testCase.verifyLessThan(min(hypot(seeds(near, 1) - c(1), seeds(near, 2) - c(2))), 1.5, 'one seed near its centre');
            testCase.verifyEqual(nnz(info.added), nnz(near) - 1, 'the rest of the grid gets no added seeds');
        end

        function testDeblendingResolvesClosePair(testCase)
            % two sources 2 px apart (less than twice their width, so their
            % sum has no dip), each varying independently across trials: the
            % dip test merges them and deblending separates them again
            src = TestFindSourcesWMF.grid();
            pair = [48.3 47.2; 48.3 49.2];
            act = TestFindSourcesWMF.simulate([src; pair], [ones(size(src, 1), 1); 2; 2], TestFindSourcesWMF.NTrials);
            ctr = mean(pair, 1);
            p = TestFindSourcesWMF.params(); p.wmfDeblendth = Inf;
            seeds = findSourcesWMF(act, p);
            testCase.verifyEqual(nnz(hypot(seeds(:, 1) - ctr(1), seeds(:, 2) - ctr(2)) < 4), 1, 'one seed without deblending');
            [seeds, info] = findSourcesWMF(act, TestFindSourcesWMF.params());
            testCase.verifyEqual(nnz(hypot(seeds(:, 1) - ctr(1), seeds(:, 2) - ctr(2)) < 4), 2, 'two seeds with deblending');
            testCase.verifyLessThan(max(TestFindSourcesWMF.nearest(pair, seeds)), 1, 'each source of the pair has a seed within 1 px');
            testCase.verifyEqual(nnz(info.added), 1, 'only the pair gets an added seed');
            testCase.verifyLessThan(max(TestFindSourcesWMF.nearest(src, seeds)), 1, 'the grid is still found');
        end

        function testFallbackKernelKeepsRowAndColumnWidths(testCase)
            % sources in pairs 5 px apart are not isolated, so K is not
            % estimated from patches; the Gaussian K it falls back to is
            % wider along columns, as these sources are
            [c, r] = meshgrid(14:20:74, 14:20:74);
            src = [r(:) c(:)] + [0.3 -0.2];
            src = [src; src + [0 5]];
            act = TestFindSourcesWMF.simulate(src, 3*ones(size(src, 1), 1), TestFindSourcesWMF.NTrials, [0.8 1.4]);
            [~, info] = findSourcesWMF(act, TestFindSourcesWMF.params());
            [X, Y] = meshgrid(-4:4);
            K = max(info.K, 0);
            testCase.verifyGreaterThan(sum(K(:).*X(:).^2)/sum(K(:).*Y(:).^2), 1.5, 'K is wider along columns than along rows');
        end

        function testUnimagedPixels(testCase)
            % columns never imaged, and a band imaged in only some trials
            src = TestFindSourcesWMF.grid();
            act = TestFindSourcesWMF.simulate(src, ones(size(src, 1), 1), TestFindSourcesWMF.NTrials);
            act(:, 1:8, :) = nan;
            act(70:78, :, 1:2:end) = nan; % half the trials: above nanThresh, left out
            seeds = findSourcesWMF(act, TestFindSourcesWMF.params());
            testCase.verifyGreaterThan(min(seeds(:, 2)), 8, 'no seed in unimaged columns');
            inBand = seeds(:, 1) >= 70 & seeds(:, 1) <= 78;
            testCase.verifyFalse(any(inBand), 'no seed where too few trials were imaged');
            keep = src(:, 2) > 12 & (src(:, 1) < 66 | src(:, 1) > 82);
            testCase.verifyLessThan(max(TestFindSourcesWMF.nearest(src(keep, :), seeds)), 1, ...
                'sources away from unimaged pixels are still found');
        end

        function testSourceActiveInOneTrialIsDropped(testCase)
            % strong enough in the combined image, but active in only one trial
            src = TestFindSourcesWMF.grid();
            act = TestFindSourcesWMF.simulate(src, ones(size(src, 1), 1), TestFindSourcesWMF.NTrials);
            c = [48 48];
            [C, R] = meshgrid(1:TestFindSourcesWMF.Sz(2), 1:TestFindSourcesWMF.Sz(1));
            act(:, :, 3) = act(:, :, 3) + 30*exp(-((R - c(1)).^2 + (C - c(2)).^2)/(2*1.2^2));
            seeds = findSourcesWMF(act, TestFindSourcesWMF.params());
            testCase.verifyGreaterThan(min(hypot(seeds(:, 1) - c(1), seeds(:, 2) - c(2))), 3);
        end

        function testOneTrialPeakDoesNotSuppressNeighbour(testCase)
            % a peak strong in the combined image but active in only one
            % trial, 4.5 px from a source active in all trials with no dip
            % between them, would merge the source into itself and then be
            % dropped; it is no candidate, so the source keeps its seed
            src = TestFindSourcesWMF.grid();
            c = [48.3 48.6]; b = c + [0 4.5];
            act = TestFindSourcesWMF.simulate([src; c], ones(size(src, 1) + 1, 1), TestFindSourcesWMF.NTrials);
            [C, R] = meshgrid(1:TestFindSourcesWMF.Sz(2), 1:TestFindSourcesWMF.Sz(1));
            act(:, :, 3) = act(:, :, 3) + 45*exp(-((R - b(1)).^2 + (C - b(2)).^2)/(2*1.2^2));
            p = TestFindSourcesWMF.params(); p.wmfDeblendth = Inf; % the merge alone, without deblending to recover the source
            seeds = findSourcesWMF(act, p);
            testCase.verifyLessThan(min(hypot(seeds(:, 1) - c(1), seeds(:, 2) - c(2))), 1, 'the source has a seed within 1 px');
            testCase.verifyGreaterThan(min(hypot(seeds(:, 1) - b(1), seeds(:, 2) - b(2))), 2, 'the one-trial peak has none');
        end

        function testExcludedPeakDoesNotSuppressNeighbour(testCase)
            % a strong source inside the exclusion mask (e.g. a soma), 4.5 px
            % from a weaker one outside it with no dip between them: the
            % masked peak is no candidate, so the weaker source keeps its
            % seed, with or without deblending, and no seed lies in the mask
            src = TestFindSourcesWMF.grid();
            s = [48.3 48.6]; a = s + [0 4.5];
            act = TestFindSourcesWMF.simulate([src; s; a], [ones(size(src, 1), 1); 3; 1], TestFindSourcesWMF.NTrials);
            [C, R] = meshgrid(1:TestFindSourcesWMF.Sz(2), 1:TestFindSourcesWMF.Sz(1));
            mask = hypot(R - s(1), C - s(2)) <= 1.5;
            p = TestFindSourcesWMF.params(); p.wmfDeblendth = Inf;
            for pp = {p, TestFindSourcesWMF.params()}
                seeds = findSourcesWMF(act, pp{1}, mask);
                testCase.verifyLessThan(min(hypot(seeds(:, 1) - a(1), seeds(:, 2) - a(2))), 1, 'the source outside the mask has a seed within 1 px');
                testCase.verifyFalse(any(mask(sub2ind(size(mask), round(seeds(:, 1)), round(seeds(:, 2))))), 'no seed in the mask');
                testCase.verifyLessThan(max(TestFindSourcesWMF.nearest(src, seeds)), 1, 'the grid is still found');
            end
        end

        function testFewTrials(testCase, nFew)
            % below 4 trials (e.g. one long recording) there are no halves to
            % test and sources need be active in only one trial; strong
            % sources, as in a long trial, are all found, and the elongated
            % source is still one seed
            src = TestFindSourcesWMF.grid();
            act = TestFindSourcesWMF.simulate(src, 4*ones(size(src, 1), 1), nFew);
            c = [48.3 48.6];
            [C, R] = meshgrid(1:TestFindSourcesWMF.Sz(2), 1:TestFindSourcesWMF.Sz(1));
            act = act + 8*TestFindSourcesWMF.Amp*exp(-((R - c(1)).^2/(2*1.1^2) + (C - c(2)).^2/(2*3^2)));
            seeds = findSourcesWMF(act, TestFindSourcesWMF.params());
            testCase.verifyLessThan(max(TestFindSourcesWMF.nearest(src, seeds)), 1, 'every source has a seed within 1 px');
            near = hypot(seeds(:, 1) - c(1), seeds(:, 2) - c(2)) < 6;
            testCase.verifyEqual(nnz(near), 1, 'one seed on the elongated source');
            testCase.verifyEqual(size(seeds, 1), size(src, 1) + 1, 'no extra seeds');
        end

        function testNoImagedPixels(testCase)
            act = nan([TestFindSourcesWMF.Sz 4]);
            [seeds, info] = findSourcesWMF(act, TestFindSourcesWMF.params());
            testCase.verifySize(seeds, [0 2]);
            testCase.verifyEmpty(info.value);
            seeds = findSourcesWMF(zeros([TestFindSourcesWMF.Sz 0]), TestFindSourcesWMF.params());
            testCase.verifySize(seeds, [0 2]);
        end
    end

    methods (Static)
        function p = params()
            p = struct('wmfPeakth', 5, 'wmfDip', 3, 'wmfDeblendth', 8, 'nanThresh', 0.33);
        end

        function src = grid()
            % sources 16 px apart, off the pixel grid, leaving the centre free
            [c, r] = meshgrid(16:16:80, 16:16:80);
            src = [r(:) c(:)] + [0.3 -0.2];
            src(hypot(src(:, 1) - 48, src(:, 2) - 48) < 1, :) = [];
        end

        function act = simulate(src, amp, nTrials, sigma)
            % rows x cols x trials activity images: row-correlated noise of
            % unit SD plus Gaussian sources (sigma 1.2 px, or [rows cols])
            % whose amplitude varies by trial
            if nargin < 4, sigma = [1.2 1.2]; end
            sz = TestFindSourcesWMF.Sz;
            [C, R] = meshgrid(1:sz(2), 1:sz(1));
            k = [1 2 1]/sqrt(6); % unit-variance smoothing along rows
            act = zeros([sz nTrials]);
            for t = 1:nTrials
                a = conv2(randn(sz + [0 2]), k, 'valid');
                for s = 1:size(src, 1)
                    g = TestFindSourcesWMF.Amp*amp(s)*(0.5 + rand);
                    a = a + g*exp(-((R - src(s, 1)).^2/(2*sigma(1)^2) + (C - src(s, 2)).^2/(2*sigma(2)^2)));
                end
                act(:, :, t) = a;
            end
        end

        function [act, c] = withElongatedSource()
            % the grid plus a source at c, 2.5x longer along columns than
            % along rows and the same in every trial
            src = TestFindSourcesWMF.grid();
            act = TestFindSourcesWMF.simulate(src, ones(size(src, 1), 1), TestFindSourcesWMF.NTrials);
            c = [48.3 48.6];
            [C, R] = meshgrid(1:TestFindSourcesWMF.Sz(2), 1:TestFindSourcesWMF.Sz(1));
            blob = exp(-((R - c(1)).^2/(2*1.1^2) + (C - c(2)).^2/(2*3^2)));
            act = act + 2*TestFindSourcesWMF.Amp*blob;
        end

        function d = nearest(src, seeds)
            % distance from each source to its nearest seed
            d = inf(size(src, 1), 1);
            for s = 1:size(src, 1)
                if ~isempty(seeds)
                    d(s) = min(hypot(seeds(:, 1) - src(s, 1), seeds(:, 2) - src(s, 2)));
                end
            end
        end
    end
end
