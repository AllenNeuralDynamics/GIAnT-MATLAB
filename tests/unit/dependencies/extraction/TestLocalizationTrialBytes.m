classdef TestLocalizationTrialBytes < matlab.unittest.TestCase
    %TESTLOCALIZATIONTRIALBYTES Unit tests for localizationTrialBytes.
    %   SILo sizes its localization pool with localizationTrialBytes. It must
    %   budget for the largest trial of a path, whichever trial that is: the
    %   first trial of a SLAP2 recording is often partial, and sizing the
    %   pool from it localized too many trials at once and ran out of memory.
    %   h5 movies are sized by their data, not the (compressed) file.
    %   Run with: runtests('tests/unit/dependencies/extraction/TestLocalizationTrialBytes.m')

    properties
        Dr
    end

    methods (TestClassSetup)
        function addPaths(testCase)
            import matlab.unittest.fixtures.PathFixture
            testsDir = fileparts(mfilename('fullpath'));
            repoRoot = fileparts(fileparts(fileparts(fileparts(testsDir))));
            testCase.applyFixture(PathFixture(fullfile(repoRoot, 'dependencies'), ...
                'IncludingSubfolders', true));
        end
    end

    methods (TestMethodSetup)
        function makeFolder(testCase)
            import matlab.unittest.fixtures.TemporaryFolderFixture
            testCase.Dr = testCase.applyFixture(TemporaryFolderFixture).Folder;
        end
    end

    methods (Test)
        function testSizesForLargestTrialNotFirst(testCase)
            % a short first trial, as at the start of a recording, then full ones
            frames = [12 300 298];
            numChannels = 2;
            fns = arrayfun(@(t) testCase.writeH5Movie(t, [30 40], frames(t), numChannels), ...
                1:numel(frames), 'UniformOutput', false);

            [bytesPerTrial, channelBytes] = localizationTrialBytes(testCase.Dr, fns, numChannels);

            testCase.verifyEqual(channelBytes, 4*30*40*frames);
            testCase.verifyEqual(bytesPerTrial, 9*4*30*40*300 + 2^30);
        end

        function testH5SizedByDataNotFile(testCase)
            % a compressed movie of constant frames is far smaller on disk
            fn = testCase.writeH5Movie(1, [64 64], 200, 1);
            dd = dir(fullfile(testCase.Dr, fn));
            testCase.assumeLessThan(dd.bytes, 4*64*64*200/10, 'the test movie should compress');

            [~, channelBytes] = localizationTrialBytes(testCase.Dr, {fn}, 1);

            testCase.verifyEqual(channelBytes, 4*64*64*200);
        end

        function testTifSizedByFile(testCase)
            fns = {'E1T1DMD1_REGISTERED_DOWNSAMPLED-80Hz.tif', 'E1T2DMD1_REGISTERED_DOWNSAMPLED-80Hz.tif'};
            nBytes = [1000 5000];
            for ix = 1:2
                fid = fopen(fullfile(testCase.Dr, fns{ix}), 'w');
                fwrite(fid, zeros(1, nBytes(ix), 'uint8'));
                fclose(fid);
            end

            [bytesPerTrial, channelBytes] = localizationTrialBytes(testCase.Dr, fns, 2);

            testCase.verifyEqual(channelBytes, nBytes/2);
            testCase.verifyEqual(bytesPerTrial, 9*5000/2 + 2^30);
        end

        function testMissingMovieErrors(testCase)
            fn = testCase.writeH5Movie(1, [8 8], 10, 1);

            testCase.verifyError(@() localizationTrialBytes(testCase.Dr, {fn, 'E1T2DMD1_missing.h5'}, 1), ...
                'localizationTrialBytes:unreadable');
        end
    end

    methods
        function fn = writeH5Movie(testCase, trial, frameSz, nFrames, numChannels)
            % as MultiRoiRegistration writes them: one /data set of frames x
            % channels, chunked by frame and compressed
            fn = sprintf('E1T%dDMD1_REGISTERED_DOWNSAMPLED-80Hz.h5', trial);
            sz = [frameSz nFrames*numChannels];
            h5create(fullfile(testCase.Dr, fn), '/data', sz, 'Datatype', 'single', ...
                'Deflate', 4, 'ChunkSize', [frameSz 1]);
            h5write(fullfile(testCase.Dr, fn), '/data', ones(sz, 'single'));
        end
    end
end
