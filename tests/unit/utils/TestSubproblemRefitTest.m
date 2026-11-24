classdef TestSubproblemRefitTest < matlab.unittest.TestCase
    %TESTSUBPROBLEMREFITTEST subproblemRefitTest on a small synthetic
    %subproblem (two sources, about one photon per pixel and sample, as on
    %SLAP2), saved as extractTrial saves one: its tables, reproducibility,
    %variants on identical data, scaled truth and injected events.
    %   Run with: runtests('tests/unit/utils/TestSubproblemRefitTest.m')

    properties
        Folder
    end

    methods (TestClassSetup)
        function addPathsAndSubproblem(testCase)
            import matlab.unittest.fixtures.PathFixture
            import matlab.unittest.fixtures.TemporaryFolderFixture
            testsDir = fileparts(mfilename('fullpath'));
            repoRoot = fileparts(fileparts(fileparts(testsDir)));
            testCase.applyFixture(PathFixture(fullfile(repoRoot, 'dependencies'), 'IncludingSubfolders', true));
            testCase.applyFixture(PathFixture(fullfile(repoRoot, 'utils')));
            testCase.Folder = testCase.applyFixture(TemporaryFolderFixture).Folder;
            TestSubproblemRefitTest.writeSubproblem(fullfile(testCase.Folder, 'path1_trial001', 'sub_001.mat'), 11);
        end
    end

    methods (Test)
        function testTables(testCase)
            R = subproblemRefitTest(testCase.Folder, 'Verbose', false);
            testCase.verifyNumElements(R.files, 1);
            testCase.verifyEqual(height(R.sources), 2);
            testCase.verifyEqual(height(R.summary), 1);
            testCase.verifyGreaterThan(height(R.events), 20);
            testCase.verifyFalse(any(R.events.injected));
            testCase.verifyTrue(all(R.events.size > 0 & R.events.trueMass > 0));
            testCase.verifyEqual(sum(R.bySize.n), height(R.events), 'every event in one size bin');
            testCase.verifyTrue(all(isfinite([R.summary.recovered R.summary.pooled R.summary.slope R.summary.F0ratio])));
            % the found events dominate: the test's regime is above the detection threshold
            testCase.verifyGreaterThan(mean(R.events.found), 0.8);
        end

        function testReproducibleAndSeeded(testCase)
            R1 = subproblemRefitTest(testCase.Folder, 'Verbose', false);
            R2 = subproblemRefitTest(testCase.Folder, 'Verbose', false);
            testCase.verifyEqual(R2.events, R1.events);
            R3 = subproblemRefitTest(testCase.Folder, 'Verbose', false, 'Seed', 5);
            testCase.verifyEqual(R3.events.trueMass, R1.events.trueMass, 'the same truth');
            testCase.verifyNotEqual(R3.events.refitMass, R1.events.refitMass, 'new data');
        end

        function testVariantsShareTheData(testCase)
            % two identical variants give identical refits; a different
            % option changes the refit only
            R = subproblemRefitTest(testCase.Folder, 'Verbose', false, ...
                'Variants', {struct(), struct(), struct('lambda', 0.3)}, 'VariantNames', {'a', 'b', 'lowLambda'});
            a = R.events(R.events.variant == "a", :);
            b = R.events(R.events.variant == "b", :);
            c = R.events(R.events.variant == "lowLambda", :);
            testCase.verifyEqual(b.refitMass, a.refitMass);
            testCase.verifyEqual(c.trueMass, a.trueMass);
            testCase.verifyNotEqual(c.refitMass, a.refitMass);
            testCase.verifyEqual(R.summary.variant, ["a"; "b"; "lowLambda"]);
        end

        function testTruthScale(testCase)
            R1 = subproblemRefitTest(testCase.Folder, 'Verbose', false);
            Rh = subproblemRefitTest(testCase.Folder, 'Verbose', false, 'TruthScale', 0.5);
            testCase.verifyEqual(Rh.events.size, 0.5 * R1.events.size, 'RelTol', 1e-12, 'the same events, half as large');
        end

        function testOverridesUpdateDerivedParams(testCase)
            % a base parameter's derived fields follow it: tau_s at its
            % saved value changes nothing, another value changes the kernel
            d = load(fullfile(testCase.Folder, 'path1_trial001', 'sub_001.mat'), 'params');
            R = subproblemRefitTest(testCase.Folder, 'Verbose', false, 'VariantNames', {'saved', 'same', 'slower'}, ...
                'Variants', {struct(), struct('tau_s', d.params.tau_s), struct('tau_s', 2*d.params.tau_s)});
            ev = @(name) R.events.refitMass(R.events.variant == name);
            testCase.verifyEqual(ev("same"), ev("saved"));
            testCase.verifyNotEqual(ev("slower"), ev("saved"));
            % tau_full alone sets the kernel too, and wins over tau_s
            R2 = subproblemRefitTest(testCase.Folder, 'Verbose', false, 'VariantNames', {'tauFull', 'both'}, ...
                'Variants', {struct('tau_full', 2*d.params.tau_full), struct('tau_s', 3*d.params.tau_s, 'tau_full', 2*d.params.tau_full)});
            testCase.verifyEqual(R2.events.refitMass(R2.events.variant == "tauFull"), ev("slower"));
            testCase.verifyEqual(R2.events.refitMass(R2.events.variant == "both"), ev("slower"));
            testCase.verifyError(@() subproblemRefitTest(testCase.Folder, 'Verbose', false, ...
                'Variants', struct('photonScale', 1)), 'subproblemRefitTest:fixedParam');
        end

        function testInjectionWithShortKernel(testCase)
            % at a low frame rate the kernel is three samples long (analyzeHz
            % 5 with tau_s 24 ms), shorter than an event's window
            R = subproblemRefitTest(testCase.Folder, 'Verbose', false, 'InjectSizes', 2, ...
                'TruthParams', struct('analyzeHz', 5), 'Variants', struct('analyzeHz', 5));
            testCase.verifyTrue(any(R.events.injected));
        end

        function testMissingFolder(testCase)
            testCase.verifyError(@() subproblemRefitTest(fullfile(testCase.Folder, 'no_such_folder'), 'Verbose', false), ...
                'subproblemRefitTest:noFiles');
        end

        function testRejectsDuplicateVariantNames(testCase)
            testCase.verifyError(@() subproblemRefitTest(testCase.Folder, 'Verbose', false, ...
                'Variants', {struct(), struct('lambda', 0.3)}, 'VariantNames', {'a', 'a'}), ...
                'subproblemRefitTest:variantNames');
        end

        function testNoEvents(testCase)
            % no true activity: empty tables with their columns, and a report
            R = subproblemRefitTest(testCase.Folder, 'Verbose', true, 'TruthScale', 0);
            testCase.verifyEqual(height(R.events), 0);
            testCase.verifyEqual(height(R.bySize), 0);
            testCase.verifyTrue(all(ismember({'variant', 'injected', 'n', 'recoveredFound'}, R.bySize.Properties.VariableNames)));
            testCase.verifyTrue(isnan(R.summary.pooled));
            % injected events alone, far below the detection threshold: none
            % is found, and recovery of found events is undefined (NaN), not 0
            R = subproblemRefitTest(testCase.Folder, 'Verbose', false, 'TruthScale', 0, 'InjectSizes', 0.02, 'InjectRate', 1);
            testCase.assertNotEmpty(R.bySize);
            testCase.verifyEqual(R.bySize.shareFound, 0);
            testCase.verifyTrue(isnan(R.bySize.recoveredFound));
        end

        function testInjectedEvents(testCase)
            sizes = [0.5 2];
            R = subproblemRefitTest(testCase.Folder, 'Verbose', false, 'InjectSizes', sizes, 'InjectRate', 1);
            inj = R.events(R.events.injected, :);
            testCase.assertNotEmpty(inj);
            testCase.verifyEqual(height(R.events(~R.events.injected, :)), ...
                height(subproblemRefitTest(testCase.Folder, 'Verbose', false).events), 'the true events are kept');
            % each injected event's size is its injected size plus what little
            % true activity is already there
            near = min(abs(inj.size ./ sizes - 1), [], 2);
            testCase.verifyLessThan(max(near), 0.5);
            testCase.verifyTrue(any(R.bySize.injected));
        end
    end

    methods (Static)
        function writeSubproblem(file, seed)
            % a subproblem as extractTrial saves one (photon-normalized data):
            % two overlapping sources, baseline about 0.8 per pixel and
            % sample, 10 s at 200 Hz, events of 1-6 units at the peak pixel
            prev = RandStream.setGlobalStream(RandStream('twister', 'Seed', seed)); % poissrnd draws from the global stream
            restore = onCleanup(@() RandStream.setGlobalStream(prev));
            params = setParams('SILo', struct('isSLAP2', true, 'nWorkers', 1));
            params.minBaseline = 0.01;
            params.numChannels = 1;
            params.selRadius = ceil(2*params.dXY);
            params.tau_full = params.tau_s*params.analyzeHz;
            params = setParamsExtractTrial(params);
            sources_p = struct('R', [7; 9], 'C', [7; 10]);
            sz = [15 16];
            selPix_p = false(sz); selPix_p(sub2ind(sz, sources_p.R, sources_p.C)) = true;
            selPix_p = imdilate(selPix_p, strel('disk', params.selRadius));
            [r, c] = find(selPix_p);
            nT = 2000;
            H = exp(-((r - sources_p.R').^2 + (c - sources_p.C').^2) / (2*1.4^2));
            H = H ./ sum(H, 1);
            sizes = 1 + 5*rand(2, nT); % at the peak pixel
            spikes = (rand(2, nT) < 0.015) .* sizes ./ max(H, [], 1)';
            X = convn(spikes, params.k / max(params.k), 'same');
            B = 0.8 + 0.2*sin((1:nT)/150) + 0.2*rand(numel(r), 1);
            F_inv_p = 0.77*ones(numel(r), nT);
            Y_p = F_inv_p .* poissrnd((B + H*X) ./ F_inv_p);
            rowSupport = [1 sz(1)]; colSupport = [1 sz(2)];
            srcIdx = [1; 2]; srcR = sources_p.R; srcC = sources_p.C; problemIx = 1; nProblems = 1; %#ok<NASGU>
            if ~isfolder(fileparts(file)), mkdir(fileparts(file)); end
            save(file, 'Y_p', 'F_inv_p', 'sources_p', 'selPix_p', 'params', 'rowSupport', 'colSupport', ...
                'srcIdx', 'srcR', 'srcC', 'problemIx', 'nProblems', '-v7.3');
        end
    end
end
