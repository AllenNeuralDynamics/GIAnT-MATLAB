classdef TestSILo < matlab.unittest.TestCase
    %TESTSILO Unit tests for SILo (non-SLAP2 path).
    %   Uses test doubles in tests/unit/source_extraction/fakes:
    %   loadAndProcessTrialAsync returns per-trial mean and activity images
    %   from a sidecar .fake.mat next to each registered file, and
    %   processAllTrials_Async records its inputs and returns minimal
    %   extraction results. Each trial's mean image is a crop of a smooth
    %   random texture at a known shift; its activity image is the same crop
    %   of Gaussian blobs plus noise, so cross-trial alignment and source
    %   positions are known. Requires Parallel Computing, Image Processing,
    %   Optimization and Statistics Toolboxes.
    %   Run with: runtests('tests/unit/source_extraction/testSILo.m')

    properties
        Dr
        SaveDr
        DataDr
        MocoDr
        Canvas
        ActCanvas
        NumTrials = 3
        NumFrames = 20
        CropSize = 40
        Margin = 5
        Shifts = [0 0; 2 -1; -1 3]          % per-trial [row col] crop shift
        Blobs = [18 22; 30 32; 34 20]       % blob centres in canvas coordinates
        ActNoise = 0.05
    end

    methods (TestClassSetup)
        function addPaths(testCase)
            import matlab.unittest.fixtures.PathFixture
            testsDir = fileparts(mfilename('fullpath'));
            repoRoot = fileparts(fileparts(fileparts(testsDir)));
            testCase.applyFixture(PathFixture(fullfile(repoRoot, 'source_extraction')));
            testCase.applyFixture(PathFixture(fullfile(repoRoot, 'dependencies'), ...
                'IncludingSubfolders', true));
            % Added last so the fakes shadow the real trial loader and extractor.
            testCase.applyFixture(PathFixture(fullfile(testsDir, 'fakes')));
        end

        function requireToolboxes(testCase)
            testCase.assumeNotEmpty(ver('images'), 'SILo requires Image Processing Toolbox.');
            testCase.assumeNotEmpty(ver('optim'), 'getActImPeaks requires Optimization Toolbox.');
            testCase.assumeNotEmpty(ver('stats'), 'getActImPeaks requires Statistics and Machine Learning Toolbox.');
            testCase.assumeNotEmpty(ver('parallel'), 'SILo requires Parallel Computing Toolbox.');

            % With nWorkers = 1 SILo deletes the pool; keep parfor from
            % starting a new one so it runs in the client with the fakes.
            ps = parallel.Settings;
            testCase.addTeardown(@restoreAutoCreate, ps, ps.Pool.AutoCreate);
            ps.Pool.AutoCreate = false;
        end
    end

    methods (TestMethodSetup)
        function makeFolders(testCase)
            import matlab.unittest.fixtures.TemporaryFolderFixture
            testCase.Dr = testCase.applyFixture(TemporaryFolderFixture).Folder;
            testCase.DataDr = testCase.applyFixture(TemporaryFolderFixture).Folder;
            testCase.SaveDr = testCase.applyFixture(TemporaryFolderFixture).Folder;
            testCase.MocoDr = fullfile(testCase.SaveDr, 'motion_correction');
            mkdir(testCase.MocoDr);
        end

        function makeCanvases(testCase)
            testCase.addTeardown(@rng, rng);
            rng(42, 'twister');
            sz = 50;
            testCase.Canvas = rescale(imgaussfilt(randn(sz), 2), 10, 100);
            [rr, cc] = ndgrid(1:sz, 1:sz);
            testCase.ActCanvas = zeros(sz);
            for b = 1:size(testCase.Blobs, 1)
                testCase.ActCanvas = testCase.ActCanvas + ...
                    exp(-((rr - testCase.Blobs(b,1)).^2 + (cc - testCase.Blobs(b,2)).^2) / 2);
            end
        end
    end

    methods (Test)
        function testAlignsTrialsToCommonFrame(testCase)
            testCase.writeDataset(1);

            testCase.runSILo(testCase.Dr);

            args = testCase.loadExtractorArgs();
            testCase.verifyEqual(args.validTrials(:)', 1:testCase.NumTrials);
            % aligned(r,c) = trial(r+m), so shift + motion is the same for all trials
            offsets = args.motOutput + testCase.Shifts';
            testCase.verifyEqual(offsets - offsets(:,1), zeros(2, testCase.NumTrials), 'AbsTol', 0.25);
        end

        function testFindsSourcesAtActivityBlobs(testCase)
            testCase.writeDataset(1);

            testCase.runSILo(testCase.Dr);

            args = testCase.loadExtractorArgs();
            src = [args.sources.R(:) args.sources.C(:)];
            testCase.assertNotEmpty(src);
            expected = testCase.expectedSourcePositions(args.motOutput);
            dist = sqrt((src(:,1) - expected(:,1)').^2 + (src(:,2) - expected(:,2)').^2);
            nBlobs = size(expected, 1);
            testCase.verifyLessThanOrEqual(min(dist, [], 1), 1.5*ones(1, nBlobs), ...
                'every blob should have a source');
            testCase.verifyLessThanOrEqual(min(dist, [], 2), 5*ones(size(src,1), 1), ...
                'sources should only be found near blobs');

            % one pixel-selection mask per source, containing the source centre
            testCase.verifyEqual(size(args.selPix, 3), size(src, 1));
            for s = 1:size(src, 1)
                testCase.verifyTrue(args.selPix(src(s,1), src(s,2), s));
            end
        end

        function testWritesSummariesAndTrialTable(testCase)
            testCase.writeDataset(1);

            testCase.runSILo(testCase.Dr);

            seDr = fullfile(testCase.SaveDr, 'source_extraction');
            testCase.verifyTrue(isfile(fullfile(seDr, 'experiment_summary.h5')));
            testCase.verifyTrue(isfile(fullfile(seDr, 'per_trial_summary.h5')));

            tt = loadStructFromH5(fullfile(testCase.Dr, 'trial_table.h5'));
            ap = tt.source_extraction.analysis_params;
            testCase.verifyEqual(ap.minBaseline, 0.01);
            testCase.verifyEqual(ap.analyzeHz, 20, 'AbsTol', 1e-9);
            testCase.verifyNotEmpty(ap.startTime);
            testCase.verifyNotEmpty(ap.endTime);
            testCase.verifyEqual(tt.source_extraction.fn_raw, ...
                arrayfun(@(t) testCase.rawName(1, t), 1:testCase.NumTrials, 'UniformOutput', false));

            args = testCase.loadExtractorArgs();
            es = loadStructFromH5(fullfile(seDr, 'experiment_summary.h5'));
            nSources = numel(args.sources.R);
            testCase.verifyEqual(double(es.Path1.visualizations.act_im_peaks), ...
                [zeros(nSources,1) args.sources.R(:)-1 args.sources.C(:)-1]);
            testCase.verifyEqual(double(es.Path1.frame_info.trial_num_frames(:)), ...
                testCase.NumFrames*ones(testCase.NumTrials, 1));
        end

        function testReturnsDerivedParams(testCase)
            testCase.writeDataset(1);

            % minBaseline is fixed by SILo and cannot be overridden
            params = testCase.runSILo(testCase.Dr, struct('minBaseline', 0.5));

            testCase.verifyEqual(params.minBaseline, 0.01);
            testCase.verifyEqual(params.numChannels, 1);
            testCase.verifyEqual(params.alignHz, 20);
            testCase.verifyEqual(params.analyzeHz, 20, 'AbsTol', 1e-9);
            testCase.verifyEqual(params.tau_full, 0.024*20, 'AbsTol', 1e-9);
            testCase.verifyEqual(params.selRadius, 6);
        end

        function testAcceptsTrialTablePath(testCase)
            testCase.writeDataset(1);
            customFn = fullfile(testCase.Dr, 'custom_table.h5');
            movefile(fullfile(testCase.Dr, 'trial_table.h5'), customFn);

            testCase.runSILo(customFn);

            tt = loadStructFromH5(customFn);
            testCase.verifyTrue(isfield(tt.source_extraction, 'analysis_params'));
            testCase.verifyFalse(isfile(fullfile(testCase.Dr, 'trial_table.h5')));
        end

        function testAcceptsJsonParams(testCase)
            testCase.writeDataset(1);
            json = '{"isSLAP2":false,"nWorkers":1,"drawUserRois":false,"photonScale":1,"peakth":10}';

            params = SILo(testCase.Dr, json);

            testCase.verifyEqual(params.peakth, 10);
            testCase.verifyEqual(params.nWorkers, 1);
            testCase.verifyFalse(params.drawUserRois);
            testCase.verifyTrue(isfile(fullfile(testCase.SaveDr, 'source_extraction', 'experiment_summary.h5')));
        end

        function testSkipsTrialWithMissingRegisteredFile(testCase)
            testCase.writeDataset(1);
            delete(fullfile(testCase.MocoDr, testCase.regName(1, 2)));

            testCase.runSILo(testCase.Dr);

            args = testCase.loadExtractorArgs();
            testCase.verifyEqual(args.validTrials(:)', [1 3]);
            testCase.verifyTrue(all(isnan(args.motOutput(:,2))));
            testCase.verifyEmpty(args.discardFrames{2});

            es = loadStructFromH5(fullfile(testCase.SaveDr, 'source_extraction', 'experiment_summary.h5'));
            testCase.verifyEqual(double(es.Path1.frame_info.trial_num_frames(:)), ...
                [testCase.NumFrames; 0; testCase.NumFrames]);
        end

        function testExcludesTrialWithFewValidActivityPixels(testCase)
            testCase.writeDataset(1);
            fakeFn = fullfile(testCase.MocoDr, [testCase.regName(1, 2) '.fake.mat']);
            fake = load(fakeFn);
            valid = false(size(fake.actIM));
            valid(1:10, 1:10) = true;
            fake.actIM(~valid) = nan;
            save(fakeFn, '-struct', 'fake');

            testCase.runSILo(testCase.Dr);

            args = testCase.loadExtractorArgs();
            testCase.verifyEqual(args.validTrials(:)', [1 3]);
            % the trial is still aligned, just not used
            testCase.verifyFalse(any(isnan(args.motOutput(:,2))));
        end

        function testSkipsImagingPathWithNoValidTrials(testCase)
            testCase.writeDataset(2);
            for t = 1:testCase.NumTrials
                delete(fullfile(testCase.MocoDr, testCase.regName(2, t)));
            end

            testCase.verifyWarning(@() testCase.runSILo(testCase.Dr), 'SILo:NoValidTrials');

            es = loadStructFromH5(fullfile(testCase.SaveDr, 'source_extraction', 'experiment_summary.h5'));
            testCase.verifyTrue(isfield(es, 'Path1'));
            testCase.verifyFalse(isfield(es, 'Path2'));
        end

        function testErrorsWhenAllTrialsMissing(testCase)
            testCase.writeDataset(1);
            for t = 1:testCase.NumTrials
                delete(fullfile(testCase.MocoDr, testCase.regName(1, t)));
            end

            try
                testCase.runSILo(testCase.Dr);
                testCase.verifyFail('Expected SILo to error when no trials have files.');
            catch ME
                testCase.verifySubstring(ME.message, 'All trials were rejected');
            end
        end

        function testNoInputErrorsWhenNonInteractive(testCase)
            testCase.assumeFalse(isInteractiveSession(), ...
                'Without inputs SILo opens a file dialog in interactive sessions.');

            testCase.verifyError(@() SILo(), 'SILo:NonInteractive');
        end
    end

    methods
        function params = runSILo(testCase, target, paramsIn)
            p = struct('isSLAP2', false, 'nWorkers', 1, 'drawUserRois', false, 'photonScale', 1);
            if nargin > 2
                for f = fieldnames(paramsIn)'
                    p.(f{1}) = paramsIn.(f{1});
                end
            end
            params = SILo(target, p);
        end

        function args = loadExtractorArgs(testCase)
            argsFn = fullfile(testCase.MocoDr, 'processAllTrials_Async_args.mat');
            testCase.assertTrue(isfile(argsFn), 'processAllTrials_Async was not called.');
            args = load(argsFn);
        end

        function pos = expectedSourcePositions(testCase, motOutput)
            % blob at canvas B is at B - Margin - shift in the trial and at
            % B - Margin - (shift + motion) in the common aligned frame
            k = round(testCase.Shifts(1,:) + motOutput(:,1)');
            pos = testCase.Blobs - testCase.Margin - k;
        end

        function writeDataset(testCase, nDMDs)
            nTrials = testCase.NumTrials;
            tt.datadr = testCase.DataDr;
            tt.savedr = testCase.SaveDr;
            tt.filename = cell(nDMDs, nTrials);
            tt.true_trial_ix = 1:nTrials;
            mc.fn_reg_ds = cell(nDMDs, nTrials);
            mc.fn_adata = cell(nDMDs, nTrials);
            mc.fn_raw = cell(nDMDs, nTrials);
            mc.registration_failed = zeros(nDMDs, nTrials);
            for d = 1:nDMDs
                for t = 1:nTrials
                    tt.filename{d,t} = sprintf('trial%d_DMD%d.tif', t, d);
                    mc.fn_reg_ds{d,t} = testCase.regName(d, t);
                    mc.fn_adata{d,t} = [testCase.stem(d, t) '_ALIGNMENTDATA.h5'];
                    mc.fn_raw{d,t} = testCase.rawName(d, t);

                    touch(fullfile(testCase.DataDr, tt.filename{d,t}));
                    touch(fullfile(testCase.MocoDr, mc.fn_reg_ds{d,t}));
                    touch(fullfile(testCase.MocoDr, mc.fn_raw{d,t}));
                    testCase.writeAlignmentData(mc.fn_adata{d,t});
                    testCase.writeTrialImages(mc.fn_reg_ds{d,t}, testCase.Shifts(t,:));
                end
            end
            tt.motion_correction = mc;
            saveStructToH5(tt, fullfile(testCase.Dr, 'trial_table.h5'));
        end

        function writeAlignmentData(testCase, fn)
            aData = struct('numChannels', 1, 'alignHz', 20, 'frametime', 0.05, ...
                'motionR', zeros(1, testCase.NumFrames), 'motionC', zeros(1, testCase.NumFrames));
            saveStructToH5(aData, fullfile(testCase.MocoDr, fn));
        end

        function writeTrialImages(testCase, regFn, shift)
            rows = testCase.Margin + (1:testCase.CropSize) + shift(1);
            cols = testCase.Margin + (1:testCase.CropSize) + shift(2);
            meanIM = testCase.Canvas(rows, cols); %#ok<NASGU>
            actIM = testCase.ActCanvas(rows, cols) + testCase.ActNoise*randn(testCase.CropSize); %#ok<NASGU>
            discardFrames = false(testCase.NumFrames, 1); %#ok<NASGU>
            peaks = struct('row', testCase.Blobs(:,1) - testCase.Margin - shift(1), ...
                'col', testCase.Blobs(:,2) - testCase.Margin - shift(2)); %#ok<NASGU>
            save(fullfile(testCase.MocoDr, [regFn '.fake.mat']), ...
                'meanIM', 'actIM', 'discardFrames', 'peaks');
        end

        function s = stem(~, d, t)
            s = sprintf('E1T%dDMD%d', t, d);
        end

        function fn = regName(testCase, d, t)
            fn = [testCase.stem(d, t) '_REGISTERED_DOWNSAMPLED-20Hz.h5'];
        end

        function fn = rawName(testCase, d, t)
            fn = [testCase.stem(d, t) '_REGISTERED_RAW.tif'];
        end
    end
end

function touch(fn)
fclose(fopen(fn, 'w'));
end

function restoreAutoCreate(ps, value)
ps.Pool.AutoCreate = value;
end
