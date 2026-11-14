classdef TestCompensatedExtraction < matlab.unittest.TestCase
    %TESTCOMPENSATEDEXTRACTION SLAP2 source extraction after registration with
    %compensateOnline. Registers a fake SLAP2 recording written as real line
    %data (a static texture under a known tissue motion, with a simulated
    %online motion correction in the line headers; see the fake
    %slap2.Slap2DataFile's lineLevel mode), then runs processAllTrials_Async's
    %SLAP2 trial loader and records the motion-corrected samples it hands to
    %source extraction (extractTrial and extractionPool are test doubles in
    %./fakes). The scene is static and noise-free, so every pixel's samples
    %should be constant; what varies is motion-correction error.
    %   Run with: runtests('tests/unit/dependencies/extraction/TestCompensatedExtraction.m')

    properties
        DataDir
        SaveDir
        TrialTableFile
        LinePeriod_s
        NumLines = 25000
        AlignHz = 80
        AnalyzeHz = 200
        MaxShift = 10
        FrameSize = [40 40]
        RoiRows = 11:30
        Name = 'acq_DMD1-TRIAL000001.dat'
    end

    methods (TestClassSetup)
        function addPaths(testCase)
            import matlab.unittest.fixtures.PathFixture
            testsDir = fileparts(mfilename('fullpath'));
            repoRoot = fileparts(fileparts(fileparts(fileparts(testsDir))));
            testCase.applyFixture(PathFixture(fullfile(repoRoot, 'motion_correction')));
            testCase.applyFixture(PathFixture(fullfile(repoRoot, 'dependencies'), 'IncludingSubfolders', true));
            testCase.applyFixture(PathFixture(fullfile(repoRoot, 'tests', 'unit', 'fakes')));
            % added last so that the doubles shadow the real extractor and pool
            testCase.applyFixture(PathFixture(fullfile(testsDir, 'fakes')));
            testCase.LinePeriod_s = slap2.Slap2DataFile.LinePeriod_s;
        end

        function requireToolboxes(testCase)
            testCase.assumeNotEmpty(ver('images'), 'Requires Image Processing Toolbox.');
            testCase.assumeNotEmpty(ver('parallel'), 'MultiRoiRegistration requires Parallel Computing Toolbox.');
        end
    end

    methods (TestMethodSetup)
        function makeFolders(testCase)
            import matlab.unittest.fixtures.TemporaryFolderFixture
            testCase.DataDir = testCase.applyFixture(TemporaryFolderFixture).Folder;
            testCase.SaveDir = testCase.applyFixture(TemporaryFolderFixture).Folder;
            testCase.TrialTableFile = fullfile(testCase.SaveDir, 'trial_table.h5');
            testCase.addTeardown(@() rmappdata0('FakeExtractTrialFile'));
        end
    end

    methods (Test)
        function testCompensatedExtractionFollowsTissue(testCase)
            % with online offset switches, the extracted samples of each pixel
            % vary much less over time when registration and extraction undo the
            % online correction
            testCase.writeAcquisition();
            [Yu, Yk] = testCase.extractBoth();
            ru = TestCompensatedExtraction.residual(Yu);
            rk = TestCompensatedExtraction.residual(Yk);
            testCase.log(1, sprintf('temporal SD / spatial SD of the extracted samples: compensated %.4f, uncompensated %.4f', rk, ru));
            testCase.verifyTrue(all(isfinite(Yk(:))));
            % part of the variation is interpolation error common to both
            % (measured: 0.034 compensated, 0.049 uncompensated)
            testCase.verifyLessThan(rk, 0.8 * ru);
            testCase.verifyLessThan(rk, 0.04);
        end

        function testExtractionWithoutSwitchesUnchanged(testCase)
            % with a constant online offset, compensation changes nothing
            testCase.writeAcquisition('updateEvery', 1e9);
            [Yu, Yk] = testCase.extractBoth();
            testCase.verifyEqual(Yk, Yu, 'AbsTol', 1e-9 * max(abs(Yu(:))));
        end
    end

    methods (Static)
        function r = residual(Y)
            % median over pixels of the temporal SD, relative to the spatial SD
            % of the pixels' mean (channel 1)
            y = double(Y(:, :, 1));
            r = median(std(y, 0, 2)) / std(mean(y, 2));
        end
    end

    methods
        function [Yu, Yk] = extractBoth(testCase)
            % register without and with compensation, then extract both at the
            % same pixels: the interior of the area both imaged
            aU = testCase.register(false);
            aK = testCase.register(true);
            testCase.verifyEqual(double(aK.slap2.trimRows), double(aU.slap2.trimRows)); % same output grid
            testCase.verifyEqual(double(aK.slap2.trimCols), double(aU.slap2.trimCols));
            valid = @(a) isfinite(squeeze(a.meanIM(1, :, :)));
            selPix = imerode(valid(aU) & valid(aK), strel('square', 7));
            Yu = testCase.extract(aU, selPix);
            Yk = testCase.extract(aK, selPix);
        end

        function aData = register(testCase, compensate)
            p = struct('nWorkers', 1, 'maxshift', testCase.MaxShift, 'alignHz', testCase.AlignHz, ...
                'saveTif', false, 'overwriteExisting', true, 'compensateOnline', compensate);
            MultiRoiRegistration(testCase.TrialTableFile, p);
            tt = loadStructFromH5(testCase.TrialTableFile);
            aData = loadStructFromH5(fullfile(testCase.SaveDir, 'motion_correction', tt.motion_correction.fn_adata{1}));
            testCase.verifyEqual(double(aData.slap2.onlineCompensation.enabled), double(compensate));
        end

        function Y = extract(testCase, aData, selPix)
            meanIM = permute(aData.meanIM, [2 3 1]); % alignment files store channels first; SILo's mean image is rows x cols x channels
            [r, c] = find(selPix);
            sources = struct('R', r, 'C', c, 'V', ones(size(r)));
            params = setParams('SILo', struct('isSLAP2', true, 'nWorkers', 1, 'analyzeHz', testCase.AnalyzeHz));
            % the parameters SILo derives before extraction
            params.minBaseline = 0.01;
            params.numChannels = double(aData.numChannels);
            params.alignHz = double(aData.alignHz);
            params.selRadius = ceil(2 * params.dXY);
            params.tau_full = params.tau_s * params.analyzeHz;
            params = setParamsExtractTrial(params);
            out = [tempname '.mat'];
            testCase.addTeardown(@() deleteIfExists(out));
            setappdata(0, 'FakeExtractTrialFile', out);
            processAllTrials_Async(testCase.DataDir, {testCase.Name}, 1, testCase.NumLines, selPix, sources, ...
                {false}, {aData}, meanIM, zeros(2, 1), {}, 1, params);
            S = load(out);
            Y = S.Y_obs;
        end

        function writeAcquisition(testCase, varargin)
            cfg = struct('frameSize', testCase.FrameSize, 'roiRows', testCase.RoiRows, ...
                'roiCols', testCase.RoiRows, 'lineLevel', true);
            for k = 1:2:numel(varargin)
                cfg.(varargin{k}) = varargin{k+1};
            end
            fn = fullfile(testCase.DataDir, testCase.Name);
            fid = fopen(fn, 'w'); fprintf(fid, '%d', testCase.NumLines); fclose(fid);
            save([fn '.fake.mat'], 'cfg');
            linePeriod_s = testCase.LinePeriod_s; %#ok<NASGU>
            save(fullfile(testCase.DataDir, 'acq_DMD1.meta'), 'linePeriod_s', '-mat');
            tt.datadr = testCase.DataDir;
            tt.savedr = testCase.SaveDir;
            tt.filename = {testCase.Name};
            tt.epoch = 1;
            tt.true_trial_ix = 1;
            tt.slap2_info.first_line = 1;
            tt.slap2_info.last_line = testCase.NumLines;
            saveStructToH5(tt, testCase.TrialTableFile);
        end
    end
end

function rmappdata0(name)
if isappdata(0, name), rmappdata(0, name); end
end

function deleteIfExists(fn)
if isfile(fn), delete(fn); end
end
