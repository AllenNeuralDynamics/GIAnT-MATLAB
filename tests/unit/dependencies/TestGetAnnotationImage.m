classdef TestGetAnnotationImage < matlab.unittest.TestCase
    %TESTGETANNOTATIONIMAGE Unit tests for getAnnotationImage's choice of trial.
    %   ROIs are annotated on the mean image of a path's longest valid trial.
    %   Movie data grows with the frames: TIFF movies are sized by their file,
    %   but h5 movies may be compressed, so their file size also depends on
    %   content and they are sized by their dataset.
    %   Run with: runtests('tests/unit/dependencies/TestGetAnnotationImage.m')

    properties
        Dr
    end

    methods (TestClassSetup)
        function addPaths(testCase)
            import matlab.unittest.fixtures.PathFixture
            testsDir = fileparts(mfilename('fullpath'));
            repoRoot = fileparts(fileparts(fileparts(testsDir)));
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
        function testH5LongestByFramesNotFileSize(testCase)
            % trial 3 is the longest, but constant frames compress far
            % better than trial 1's noise
            frames = [100 100 300];
            noisy = [true false false];
            tt = testCase.writeH5Trials(frames, noisy, true);
            testCase.assumeGreaterThan(testCase.fileBytes(tt, 1), testCase.fileBytes(tt, 3), ...
                'the noisy movie should be the largest file');

            [IM, fn] = getAnnotationImage(tt, true(1, 3), testCase.Dr, 1);

            testCase.verifyEqual(fn, 'E1T3DMD1_REGISTERED_DOWNSAMPLED-80Hz');
            testCase.verifyEqual(IM, 3*ones(30, 40, 'single'));
        end

        function testH5TieGoesToFirstTrial(testCase)
            % equal trials of a continuous recording choose the same trial
            % as their TIFF movies would
            tt = testCase.writeH5Trials([100 100 100], [false true false], true);
            testCase.assumeGreaterThan(testCase.fileBytes(tt, 2), testCase.fileBytes(tt, 1), ...
                'the noisy movie should be the largest file');

            [~, fn] = getAnnotationImage(tt, true(1, 3), testCase.Dr, 1);

            testCase.verifyEqual(fn, 'E1T1DMD1_REGISTERED_DOWNSAMPLED-80Hz');
        end

        function testKeepTrialsAndMissingMovies(testCase)
            tt = testCase.writeH5Trials([100 300 200 100], false(1, 4), true);
            delete(fullfile(testCase.Dr, tt.motion_correction.fn_reg_ds{1, 3}));

            % trial 2 is excluded and trial 3's movie is missing
            [~, fn] = getAnnotationImage(tt, logical([1 0 1 1]), testCase.Dr, 1);

            testCase.verifyEqual(fn, 'E1T1DMD1_REGISTERED_DOWNSAMPLED-80Hz');
        end

        function testTifSizedByFile(testCase)
            nBytes = [1000 5000 5000];
            tt = testCase.trialTable(3, '.tif');
            for t = 1:3
                fid = fopen(fullfile(testCase.Dr, tt.motion_correction.fn_reg_ds{1, t}), 'w');
                fwrite(fid, zeros(1, nBytes(t), 'uint8'));
                fclose(fid);
                testCase.writeAlignment(tt, t, true);
            end

            [IM, fn] = getAnnotationImage(tt, true(1, 3), testCase.Dr, 1);

            testCase.verifyEqual(fn, 'E1T2DMD1_REGISTERED_DOWNSAMPLED-80Hz');
            testCase.verifyEqual(IM, 2*ones(30, 40, 'single'));
        end

        function testH5MovieWithoutMeanIM(testCase)
            tt = testCase.writeH5Trials([10 20], [false true], false);

            [IM, fn] = getAnnotationImage(tt, true(1, 2), testCase.Dr, 1);

            testCase.verifyEqual(fn, 'E1T2DMD1_REGISTERED_DOWNSAMPLED-80Hz');
            data = h5read(fullfile(testCase.Dr, [fn '.h5']), '/data');
            testCase.verifyEqual(IM, mean(data, 3, 'omitnan'), 'AbsTol', 1e-6);
        end

        function testNoValidTrials(testCase)
            tt = testCase.trialTable(2, '.h5');

            [IM, fn] = getAnnotationImage(tt, false(1, 2), testCase.Dr, 1);

            testCase.verifyEmpty(IM);
            testCase.verifyEmpty(fn);
        end
    end

    methods
        function tt = trialTable(~, nTrials, ext)
            stems = arrayfun(@(t) sprintf('E1T%dDMD1', t), 1:nTrials, 'UniformOutput', false);
            tt.motion_correction.fn_reg_ds = strcat(stems, ['_REGISTERED_DOWNSAMPLED-80Hz' ext]);
            tt.motion_correction.fn_adata = strcat(stems, '_ALIGNMENTDATA.h5');
        end

        function tt = writeH5Trials(testCase, frames, noisy, withMeanIM)
            % movies as MultiRoiRegistration writes them with saveTif=false:
            % one /data set of 2-channel frames, chunked by frame, compressed
            numChannels = 2;
            frameSz = [30 40];
            tt = testCase.trialTable(numel(frames), '.h5');
            rng(0);
            for t = 1:numel(frames)
                fn = fullfile(testCase.Dr, tt.motion_correction.fn_reg_ds{1, t});
                sz = [frameSz frames(t)*numChannels];
                h5create(fn, '/data', sz, 'Datatype', 'single', 'Deflate', 4, ...
                    'ChunkSize', [frameSz 1]);
                if noisy(t)
                    h5write(fn, '/data', rand(sz, 'single'));
                else
                    h5write(fn, '/data', ones(sz, 'single'));
                end
                testCase.writeAlignment(tt, t, withMeanIM);
            end
        end

        function writeAlignment(testCase, tt, t, withMeanIM)
            % meanIM is numChannels x H x W; its value identifies the trial
            s.alignHz = 80;
            if withMeanIM
                s.meanIM = t*ones(2, 30, 40, 'single');
            end
            saveStructToH5(s, fullfile(testCase.Dr, tt.motion_correction.fn_adata{1, t}));
        end

        function n = fileBytes(testCase, tt, t)
            d = dir(fullfile(testCase.Dr, tt.motion_correction.fn_reg_ds{1, t}));
            n = d.bytes;
        end
    end
end
