classdef TestLoadAndProcessTrialAsync < matlab.unittest.TestCase
    %TESTLOADANDPROCESSTRIALASYNC Unit tests for loadAndProcessTrialAsync.
    %   loadAndProcessTrialAsync loads a trial's registered movie and hands
    %   the activity channel and its variance factors to localizeSources_vIM
    %   (replaced here by the recording fake in localizationFakes). It keeps
    %   no copy of the movie while localizing, so if localization fails it
    %   loads the movie again and retries with the image edges excluded.
    %   Run with: runtests('tests/unit/dependencies/extraction/TestLoadAndProcessTrialAsync.m')

    properties
        Dr
        Movie       % rows x cols x channels x frames
        VarFactors  % rows x cols x frames
        Params
        Fn = 'E1T1DMD1_REGISTERED_DOWNSAMPLED-20Hz.h5'
        AlignHz = 20
    end

    methods (TestClassSetup)
        function addPaths(testCase)
            import matlab.unittest.fixtures.PathFixture
            testsDir = fileparts(mfilename('fullpath'));
            repoRoot = fileparts(fileparts(fileparts(fileparts(testsDir))));
            testCase.applyFixture(PathFixture(fullfile(repoRoot, 'dependencies'), ...
                'IncludingSubfolders', true));
            testCase.applyFixture(PathFixture(fullfile(testsDir, 'localizationFakes')));
        end
    end

    methods (TestMethodSetup)
        function writeTrial(testCase)
            import matlab.unittest.fixtures.TemporaryFolderFixture
            testCase.Dr = testCase.applyFixture(TemporaryFolderFixture).Folder;
            rng(3);
            rows = 12; cols = 10; numChannels = 2; nFrames = 40;
            testCase.Movie = rand(rows, cols, numChannels, nFrames, 'single');
            testCase.Movie(:, 1:2, :, :) = nan; % not imaged
            testCase.VarFactors = 1 + rand(rows, cols, nFrames, 'single');

            % interleaved frames x channels, as MultiRoiRegistration writes them
            sz = [rows cols numChannels*nFrames];
            h5create(fullfile(testCase.Dr, testCase.Fn), '/data', sz, 'Datatype', 'single');
            h5write(fullfile(testCase.Dr, testCase.Fn), '/data', reshape(testCase.Movie, sz));
            aData = struct('alignHz', testCase.AlignHz, 'recNegErr', zeros(1, nFrames), ...
                'slap2', struct('varFacDS', testCase.VarFactors));
            saveStructToH5(aData, fullfile(testCase.Dr, 'E1T1DMD1_ALIGNMENTDATA.h5'));

            testCase.Params = struct('nanThresh', 0.5, 'activityChannel', 2, ...
                'discardInitial_s', 0.1, 'motionThresh', 3, 'interTrialAlignment', 'rigid');
            setappdata(groot, 'localizeCalls', {});
            setappdata(groot, 'localizeFailures', 0);
            testCase.addTeardown(@rmappdata, groot, 'localizeCalls');
            testCase.addTeardown(@rmappdata, groot, 'localizeFailures');
        end
    end

    methods (Test)
        function testLocalizesActivityChannel(testCase)
            [meanIM, IMc, aData, peaks, discardFrames] = loadAndProcessTrialAsync( ...
                testCase.Dr, testCase.Fn, 2, testCase.Params);

            calls = getappdata(groot, 'localizeCalls');
            testCase.verifyNumElements(calls, 1);
            testCase.verifyEqual(calls{1}.IM, testCase.activityMovie(discardFrames));
            testCase.verifyEqual(calls{1}.vIM, testCase.VarFactors);
            testCase.verifyEqual(calls{1}.alignHz, testCase.AlignHz);
            testCase.verifyEqual(find(discardFrames(:))', 1:2, 'the first 0.1 s is discarded');
            testCase.verifyFalse(isfield(aData.slap2, 'varFacDS'), ...
                'the variance factors are handed over, not returned');
            testCase.verifyEqual(meanIM, mean(testCase.Movie, 4, 'omitnan'));
            testCase.verifyEqual(IMc, mean(calls{1}.IM, 3, 'omitnan'));
            testCase.verifyEqual(peaks.row, 2);
        end

        function testRetriesWithEdgesExcluded(testCase)
            setappdata(groot, 'localizeFailures', 1);

            [~, IMc, ~, peaks, discardFrames] = loadAndProcessTrialAsync( ...
                testCase.Dr, testCase.Fn, 2, testCase.Params);

            calls = getappdata(groot, 'localizeCalls');
            testCase.verifyNumElements(calls, 2);
            expected = testCase.activityMovie(discardFrames);
            testCase.verifyEqual(calls{1}.IM, expected);
            expected([1 end], :, :) = nan;
            expected(:, [1 end], :) = nan;
            testCase.verifyEqual(calls{2}.IM, expected, ...
                'the retry gets the movie, reloaded, with the edges excluded');
            testCase.verifyEqual(calls{2}.vIM, testCase.VarFactors);
            testCase.verifyEqual(IMc, mean(expected, 3, 'omitnan'));
            testCase.verifyEqual(peaks.row, 2);
        end

        function testRetryFailureIsRaised(testCase)
            setappdata(groot, 'localizeFailures', 2);

            testCase.verifyError(@() loadAndProcessTrialAsync(testCase.Dr, testCase.Fn, 2, testCase.Params), ...
                'localizeSources_vIM:fake');
        end
    end

    methods
        function IM = activityMovie(testCase, discardFrames)
            IM = squeeze(testCase.Movie(:, :, testCase.Params.activityChannel, :));
            IM(:, :, discardFrames) = nan;
        end
    end
end
