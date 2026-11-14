classdef TestMultiRoiRegistration < matlab.unittest.TestCase
    %TESTMULTIROIREGISTRATION Unit tests for MultiRoiRegistration.
    %   Uses test doubles in tests/unit/fakes for the SLAP2 reader and the
    %   BigTIFF writer. Each fake recording is a smooth random texture in a
    %   fixed ROI, optionally shifted part-way through to simulate motion.
    %   Requires Parallel Computing Toolbox, because MultiRoiRegistration
    %   always queries (and may start) a process pool.
    %   Run with: runtests('tests/unit/motion_correction/TestMultiRoiRegistration.m')

    properties
        DataDir
        SaveDir
        TrialTableFile
        LinePeriod_s
        NumLines = 25000
        AlignHz = 80
        MaxShift = 10
        FrameSize = [40 40]
        RoiRows = 11:30
    end

    methods (TestClassSetup)
        function addPaths(testCase)
            import matlab.unittest.fixtures.PathFixture
            testsDir = fileparts(mfilename('fullpath'));
            unitDir = fileparts(testsDir);
            repoRoot = fileparts(fileparts(unitDir));
            testCase.applyFixture(PathFixture(fullfile(repoRoot, 'motion_correction')));
            testCase.applyFixture(PathFixture(fullfile(repoRoot, 'dependencies'), ...
                'IncludingSubfolders', true));
            % Added last so the fakes shadow the real reader and writer.
            testCase.applyFixture(PathFixture(fullfile(unitDir, 'fakes')));
            testCase.LinePeriod_s = slap2.Slap2DataFile.LinePeriod_s;
        end

        function requireParallelToolbox(testCase)
            testCase.assumeNotEmpty(ver('parallel'), ...
                'MultiRoiRegistration requires Parallel Computing Toolbox.');
            if isempty(gcp('nocreate'))
                % Shut down the pool MultiRoiRegistration starts.
                testCase.addTeardown(@() delete(gcp('nocreate')));
            end
        end
    end

    methods (TestMethodSetup)
        function makeFolders(testCase)
            import matlab.unittest.fixtures.TemporaryFolderFixture
            testCase.DataDir = testCase.applyFixture(TemporaryFolderFixture).Folder;
            testCase.SaveDir = testCase.applyFixture(TemporaryFolderFixture).Folder;
            testCase.TrialTableFile = fullfile(testCase.SaveDir, 'trial_table.h5');
        end
    end

    methods (Test)
        function testWriteAlignmentFiles(testCase)
            names = {'acq_DMD1-TRIAL000001.dat', 'acq_DMD1-TRIAL000002.dat'; ...
                     'acq_DMD2-TRIAL000001.dat', 'acq_DMD2-TRIAL000002.dat'};
            for k = 1:numel(names)
                testCase.writeAcquisition(names{k});
            end
            testCase.writeTrialTable(names);

            [tt, params] = testCase.runRegistration();

            mc = tt.motion_correction;

            % verify the correct files were output
            [regNames, adataNames] = testCase.outputNames(2, 2, '.h5');
            testCase.verifyEqual(mc.fn_reg_ds, regNames);
            testCase.verifyEqual(mc.fn_adata, adataNames);
            testCase.verifyEqual(double(mc.registration_failed), zeros(2, 2));
            for k = 1:numel(regNames)
                testCase.verifyTrue(isfile(testCase.mocoFile(regNames{k})), regNames{k});
                testCase.verifyTrue(isfile(testCase.mocoFile(adataNames{k})), adataNames{k});
            end

            testCase.verifyEqual(params.maxshift, testCase.MaxShift);
            testCase.verifyEqual(params.clipShift, 5); % default kept
            testCase.verifyNotEmpty(params.startTime);
            testCase.verifyNotEmpty(params.endTime);
            testCase.verifyEqual(mc.align_params.maxshift, testCase.MaxShift);

            % verify the contents of the alignment file
            aData = testCase.loadAlignmentData(tt);
            expectedFrames = testCase.expectedDSframes();
            nFrames = numel(expectedFrames);
            padded = numel(testCase.RoiRows) + 2*testCase.MaxShift;
            testCase.verifyEqual(double(aData.DSframes(:)'), expectedFrames);
            testCase.verifyEqual(aData.frametime, 1/testCase.AlignHz);
            testCase.verifyEqual(double(aData.numChannels), 1);
            testCase.verifyEqual(aData.motionDSr, zeros(1, nFrames), 'AbsTol', 0.25);
            testCase.verifyEqual(aData.motionDSc, zeros(1, nFrames), 'AbsTol', 0.25);
            testCase.verifyEqual(double(aData.registrationFailed), 0);
            testCase.verifySize(aData.meanIM, [1 padded padded]);
            testCase.verifyEqual(double(aData.slap2.trimRows(:)'), testCase.RoiRows);
            testCase.verifyEqual(double(aData.slap2.trimCols(:)'), testCase.RoiRows);
            testCase.verifyEqual(double(aData.slap2.cropRow), testCase.RoiRows(1) - testCase.MaxShift);
            testCase.verifyEqual(double(aData.slap2.cropCol), testCase.RoiRows(1) - testCase.MaxShift);
        end

        function testRecNegErrForShortRecording(testCase)
            % A recording shorter than one 10 s chunk still gets a
            % reconstruction error for every frame (it used to get none).
            name = 'acq_DMD1-TRIAL000001.dat';
            testCase.writeAcquisition(name);
            testCase.writeTrialTable({name});

            tt = testCase.runRegistration();

            aData = testCase.loadAlignmentData(tt);
            nFrames = numel(testCase.expectedDSframes());
            testCase.assertLessThan(nFrames, 10*testCase.AlignHz);
            testCase.verifyNumElements(aData.recNegErr, nFrames);
            testCase.verifyTrue(all(isfinite(aData.recNegErr(:))));
        end

        function testRecoverKnownShift(testCase)
            name = 'acq_DMD1-TRIAL000001.dat';
            shiftLine = 0.75*testCase.NumLines;
            testCase.writeAcquisition(name, 'shiftAfterLine', shiftLine, 'shift', [2 -1]);
            testCase.writeTrialTable({name});

            tt = testCase.runRegistration();

            aData = testCase.loadAlignmentData(tt);
            after = aData.DSframes >= shiftLine;
            testCase.verifyEqual(median(aData.motionDSr(~after)), 0, 'AbsTol', 0.25);
            testCase.verifyEqual(median(aData.motionDSc(~after)), 0, 'AbsTol', 0.25);
            % The fake samples the texture at +shift, so registration must
            % move the frame by -shift to line it up again.
            testCase.verifyEqual(median(aData.motionDSr(after)), -2, 'AbsTol', 0.25);
            testCase.verifyEqual(median(aData.motionDSc(after)), 1, 'AbsTol', 0.25);
            testCase.verifyEqual(double(aData.registrationFailed), 0);
        end

        function testTwoChannelsInterleaved(testCase)
            name = 'acq_DMD1-TRIAL000001.dat';
            testCase.writeAcquisition(name, 'numChannels', 2);
            testCase.writeTrialTable({name});

            tt = testCase.runRegistration();

            aData = testCase.loadAlignmentData(tt);
            nFrames = numel(testCase.expectedDSframes());
            padded = numel(testCase.RoiRows) + 2*testCase.MaxShift;
            testCase.verifyEqual(double(aData.numChannels), 2);
            testCase.verifySize(aData.meanIM, [2 padded padded]);

            regFile = testCase.mocoFile(tt.motion_correction.fn_reg_ds{1});
            info = h5info(regFile, '/data');
            testCase.verifyEqual(info.Dataspace.Size, [padded padded 2*nFrames]);
            % The fake's channel 2 is channel 1 at half intensity, so pages
            % should alternate ch1, ch2 within each frame.
            pages = h5read(regFile, '/data', [1 1 1], [padded padded 2]);
            testCase.verifyEqual(pages(:,:,2), pages(:,:,1)/2, 'RelTol', 1e-5);
        end

        function testSaveTif(testCase)
            name = 'acq_DMD1-TRIAL000001.dat';
            testCase.writeAcquisition(name);
            testCase.writeTrialTable({name});

            tt = testCase.runRegistration('saveTif', true);

            regNames = testCase.outputNames(1, 1, '.tif');
            testCase.verifyEqual(tt.motion_correction.fn_reg_ds, regNames);
            % The fake TIFF writer stores raw single pages (4 bytes per pixel).
            padded = numel(testCase.RoiRows) + 2*testCase.MaxShift;
            nFrames = numel(testCase.expectedDSframes());
            d = dir(testCase.mocoFile(regNames{1}));
            testCase.verifyEqual(d.bytes, padded*padded*4*nFrames);
        end

        function testSkipExistingOutputs(testCase)
            % No recording on disk, so any alignment attempt would error.
            name = 'acq_DMD1-TRIAL000001.dat';
            testCase.writeTrialTable({name});
            [regNames, adataNames] = testCase.outputNames(1, 1, '.h5');
            mkdir(fullfile(testCase.SaveDir, 'motion_correction'));
            testCase.writeText(testCase.mocoFile(regNames{1}), 'existing');
            testCase.writeText(testCase.mocoFile(adataNames{1}), 'existing');

            tt = testCase.runRegistration('overwriteExisting', false);

            testCase.verifyEqual(fileread(testCase.mocoFile(regNames{1})), 'existing');
            testCase.verifyEqual(fileread(testCase.mocoFile(adataNames{1})), 'existing');
            testCase.verifyEqual(tt.motion_correction.fn_reg_ds, regNames);
            testCase.verifyEqual(tt.motion_correction.fn_adata, adataNames);
            testCase.verifyEqual(double(tt.motion_correction.registration_failed), 0);
        end

        function testOverwriteExistingOutputs(testCase)
            name = 'acq_DMD1-TRIAL000001.dat';
            testCase.writeAcquisition(name);
            testCase.writeTrialTable({name});
            [regNames, adataNames] = testCase.outputNames(1, 1, '.h5');
            mkdir(fullfile(testCase.SaveDir, 'motion_correction'));
            testCase.writeText(testCase.mocoFile(regNames{1}), 'existing');
            testCase.writeText(testCase.mocoFile(adataNames{1}), 'existing');

            tt = testCase.runRegistration('overwriteExisting', true);

            aData = testCase.loadAlignmentData(tt);
            testCase.verifyTrue(isfield(aData, 'motionDSr'));
            info = h5info(testCase.mocoFile(regNames{1}), '/data');
            testCase.verifyEqual(info.Dataspace.Size(3), numel(testCase.expectedDSframes()));
        end

        function testReadErrorFlagsTrial(testCase)
            % 40 template frames are read first; fail on the first
            % registration frame, inside MultiRoiRegistration's try/catch.
            name = 'acq_DMD1-TRIAL000001.dat';
            testCase.writeAcquisition(name, 'failOnCall', 41);
            testCase.writeTrialTable({name});

            tt = testCase.runRegistration();

            testCase.verifyEqual(double(tt.motion_correction.registration_failed), 1);
            testCase.verifyFalse(isfile(testCase.mocoFile(tt.motion_correction.fn_adata{1})));
        end

        function testReVoltUsesFirstChannel(testCase)
            name = 'acq_DMD1-TRIAL000001.dat';
            testCase.writeAcquisition(name, 'numChannels', 2);
            % first_line_original marks first_line as already adjusted, so
            % laser-onset detection is skipped.
            testCase.writeTrialTable({name}, struct('first_line_original', 1));

            tt = testCase.runRegistration('isReVolt', true);

            aData = testCase.loadAlignmentData(tt);
            testCase.verifyEqual(double(aData.numChannels), 1);
            info = h5info(testCase.mocoFile(tt.motion_correction.fn_reg_ds{1}), '/data');
            testCase.verifyEqual(info.Dataspace.Size(3), numel(testCase.expectedDSframes()));
            testCase.verifyEqual(double(tt.slap2_info.first_line), 1);
            testCase.verifyEqual(double(tt.motion_correction.first_line_original), 1);
        end

        function testMaxShiftExceededError(testCase)
            name = 'acq_DMD1-TRIAL000001.dat';
            testCase.writeAcquisition(name, ...
                'shiftAfterLine', 0.5*testCase.NumLines, 'shift', [4 0]);
            testCase.writeTrialTable({name});

            tt = testCase.runRegistration();

            % Flagged, but the alignment data is still saved.
            aData = testCase.loadAlignmentData(tt);
            testCase.verifyGreaterThan(std(aData.motionDSr), 1.5);
            testCase.verifyEqual(double(aData.registrationFailed), 1);
            testCase.verifyEqual(double(tt.motion_correction.registration_failed), 1);
        end

        function testReVoltRefStackTemplateError(testCase)
            name = 'acq_DMD1-TRIAL000001.dat';
            testCase.writeAcquisition(name);
            testCase.writeTrialTable({name}, struct('first_line_original', 1));

            testCase.verifyError(@() testCase.runRegistration( ...
                'isReVolt', true, 'refStackTemplate', true), ?MException);
        end

        function testNonIntegerMaxshiftError(testCase)
            name = 'acq_DMD1-TRIAL000001.dat';
            testCase.writeAcquisition(name);
            testCase.writeTrialTable({name});

            testCase.verifyError(@() testCase.runRegistration('maxshift', 2.5), ...
                ?MException);
        end
    end

    methods (Access = private)
        function [tt, params] = runRegistration(testCase, varargin)
            paramsIn = struct('nWorkers', 1, 'maxshift', testCase.MaxShift, ...
                'alignHz', testCase.AlignHz, 'saveTif', false);
            for k = 1:2:numel(varargin)
                paramsIn.(varargin{k}) = varargin{k+1};
            end
            params = MultiRoiRegistration(testCase.TrialTableFile, paramsIn);
            tt = loadStructFromH5(testCase.TrialTableFile);
        end

        function writeAcquisition(testCase, name, varargin)
            % Writes a fake recording: the .dat line count, the image config
            % sidecar read by the fake reader, and the per-DMD .meta file.
            cfg = struct('frameSize', testCase.FrameSize, ...
                'roiRows', testCase.RoiRows, 'roiCols', testCase.RoiRows);
            for k = 1:2:numel(varargin)
                cfg.(varargin{k}) = varargin{k+1};
            end
            fn = fullfile(testCase.DataDir, name);
            testCase.writeText(fn, sprintf('%d', testCase.NumLines));
            save([fn '.fake.mat'], 'cfg');
            linePeriod_s = testCase.LinePeriod_s;
            metaName = [regexp(name, '^.*?DMD\d', 'match', 'once') '.meta'];
            save(fullfile(testCase.DataDir, metaName), 'linePeriod_s', '-mat');
        end

        function writeTrialTable(testCase, names, motionCorrection)
            sz = size(names);
            tt.datadr = testCase.DataDir;
            tt.savedr = testCase.SaveDir;
            tt.filename = names;
            tt.epoch = ones(sz);
            tt.true_trial_ix = repmat(1:sz(2), sz(1), 1);
            tt.slap2_info.first_line = ones(sz);
            tt.slap2_info.last_line = testCase.NumLines*ones(sz);
            if nargin > 2
                tt.motion_correction = motionCorrection;
            end
            saveStructToH5(tt, testCase.TrialTableFile);
        end

        function [regNames, adataNames] = outputNames(testCase, nDMDs, nTrials, ext)
            regNames = cell(nDMDs, nTrials);
            adataNames = cell(nDMDs, nTrials);
            for d = 1:nDMDs
                for t = 1:nTrials
                    regNames{d, t} = sprintf('E1T%dDMD%d_REGISTERED_DOWNSAMPLED-%dHz%s', ...
                        t, d, testCase.AlignHz, ext);
                    adataNames{d, t} = sprintf('E1T%dDMD%d_ALIGNMENTDATA.h5', t, d);
                end
            end
        end

        function frames = expectedDSframes(testCase)
            % Mirrors how MultiRoiRegistration picks downsampled frame times.
            dt = (1/testCase.LinePeriod_s)/testCase.AlignHz;
            frames = ceil(1:dt:testCase.NumLines);
        end

        function aData = loadAlignmentData(testCase, tt)
            aData = loadStructFromH5(testCase.mocoFile(tt.motion_correction.fn_adata{1}));
        end

        function fn = mocoFile(testCase, name)
            fn = fullfile(testCase.SaveDir, 'motion_correction', name);
        end

        function writeText(~, fn, contents)
            fid = fopen(fn, 'w');
            fprintf(fid, '%s', contents);
            fclose(fid);
        end
    end
end
