classdef TestBuildTrialTableSLAP2 < matlab.unittest.TestCase
    %TESTBUILDTRIALTABLESLAP2 Unit tests for buildTrialTableSLAP2.
    %   Uses test doubles in tests/unit/fakes for the SLAP2 reader and the
    %   reference-stack TIFF reader, so no real SLAP2 data is needed.
    %   Run with: runtests('tests/unit/TestBuildTrialTableSLAP2.m')

    properties
        DataDir
        SaveDir
        LinePeriod_s
        StartTime = datetime(2024, 2, 1, 12, 0, 0, 'TimeZone', 'local')
    end

    methods (TestClassSetup)
        function addPaths(testCase)
            import matlab.unittest.fixtures.PathFixture
            testsDir = fileparts(mfilename('fullpath'));
            repoRoot = fileparts(fileparts(testsDir));
            testCase.applyFixture(PathFixture(repoRoot));
            testCase.applyFixture(PathFixture(fullfile(repoRoot, 'dependencies'), ...
                'IncludingSubfolders', true));
            % Added last so the fakes shadow the real readers.
            testCase.applyFixture(PathFixture(fullfile(testsDir, 'fakes')));
            testCase.LinePeriod_s = slap2.Slap2DataFile.LinePeriod_s;
        end
    end

    methods (TestMethodSetup)
        function makeFolders(testCase)
            import matlab.unittest.fixtures.TemporaryFolderFixture
            testCase.DataDir = testCase.applyFixture(TemporaryFolderFixture).Folder;
            testCase.SaveDir = testCase.applyFixture(TemporaryFolderFixture).Folder;
        end
    end

    methods (Test)
        function testTriggeredTrials(testCase)
            nLines = [100000 120000 110000];
            [dmd1, dmd2] = testCase.writeTriggeredFiles(nLines, nLines);
            testCase.writeReferenceStacks();

            tt = buildTrialTableSLAP2(testCase.DataDir, testCase.SaveDir, true);

            testCase.verifyEqual(tt.filename, [dmd1; dmd2]);
            testCase.verifyEqual(tt.epoch, ones(2, 3));
            testCase.verifyEqual(tt.true_trial_ix, repmat(1:3, 2, 1));
            testCase.verifyEqual(tt.slap2_info.first_line, ones(2, 3));
            testCase.verifyEqual(tt.slap2_info.last_line, repmat(nLines, 2, 1));

            endTimes = testCase.datenums(dmd1);
            testCase.verifyEqual(tt.slap2_info.trial_end_time_from_pc, endTimes, ...
                'AbsTol', 1e-9);
            testCase.verifyEqual(tt.slap2_info.trial_start_time_inferred, ...
                endTimes - nLines*testCase.LinePeriod_s/86400, 'AbsTol', 1e-9);
        end

        function testTriggeredTrialLineMismatch(testCase)
            % DMD1 failed to increment its trial number, so its second file
            % holds DMD2's second and third trials back to back.
            [dmd1, dmd2] = testCase.writeTriggeredFiles( ...
                [100000 220000], [100000 100000 120000]);
            testCase.writeReferenceStacks();

            tt = buildTrialTableSLAP2(testCase.DataDir, testCase.SaveDir, true);

            testCase.verifyEqual(tt.filename, ...
                [dmd1([1 2 2]); dmd2([1 2 3])]);
            testCase.verifyEqual(tt.slap2_info.first_line, ...
                [1 1 100001; 1 1 1]);
            testCase.verifyEqual(tt.slap2_info.last_line, ...
                [100000 100000 220000; 100000 100000 120000]);

            dn1 = testCase.datenums(dmd1);
            dn2 = testCase.datenums(dmd2);
            testCase.verifyEqual(tt.slap2_info.trial_end_time_from_pc, ...
                [dn1(1) dn2(2) dn1(2)], 'AbsTol', 1e-9);
        end

        function testIgnoreMulticycleFiles(testCase)
            nLines = 450000;
            names = {'acq_DMD1-CYCLE-000000.dat', 'acq_DMD2-CYCLE-000000.dat'};
            testCase.writeDat(names{1}, nLines, testCase.StartTime);
            testCase.writeDat(names{2}, nLines, testCase.StartTime);
            % Later cycles must be skipped; unparseable contents would
            % produce NaN line counts if they were read.
            testCase.writeText('acq_DMD1-CYCLE-000001.dat', 'not a line count');
            testCase.writeText('acq_DMD2-CYCLE-000001.dat', 'not a line count');
            testCase.writeReferenceStacks();

            tt = buildTrialTableSLAP2(testCase.DataDir, testCase.SaveDir, true);

            edges = [1 150001 300001 450001];
            testCase.verifyEqual(tt.filename, repmat(names', 1, 3));
            testCase.verifyEqual(tt.true_trial_ix, repmat(1:3, 2, 1));
            testCase.verifyEqual(tt.slap2_info.first_line, repmat(edges(1:3), 2, 1));
            testCase.verifyEqual(tt.slap2_info.last_line, repmat(edges(2:4), 2, 1));

            dn = testCase.datenums(names(1));
            toDays = testCase.LinePeriod_s/86400;
            testCase.verifyEqual(tt.slap2_info.trial_end_time_from_pc, ...
                dn - (nLines - edges(2:4))*toDays, 'AbsTol', 1e-9);
            testCase.verifyEqual(tt.slap2_info.trial_start_time_inferred, ...
                dn - (nLines - edges(1:3))*toDays, 'AbsTol', 1e-9);
        end

        function testWriteReferenceStack(testCase)
            testCase.writeTriggeredFiles(100000, 100000);
            testCase.writeReferenceStacks();

            tt = buildTrialTableSLAP2(testCase.DataDir, testCase.SaveDir, true);

            for p = {'Path1', 'Path2'}
                ref = tt.slap2_info.ref_stack.(p{1});
                testCase.verifyEqual(ref.channels, [1 2]);
                testCase.verifyEqual(ref.Zs, [0 1 2]);
                testCase.verifyEqual(ref.dmdPixel2SampleTransform, eye(3));
                testCase.verifySize(ref.IM, [8 8 6]);
            end
        end

        function testWriteTrialTableH5(testCase)
            [dmd1, dmd2] = testCase.writeTriggeredFiles([100000 120000], [100000 120000]);
            testCase.writeReferenceStacks();

            buildTrialTableSLAP2(testCase.DataDir, testCase.SaveDir, true);

            h5fn = fullfile(testCase.SaveDir, 'trial_table.h5');
            testCase.assertTrue(isfile(h5fn), 'trial_table.h5 was not written');
            tt = loadStructFromH5(h5fn);
            testCase.verifyEqual(tt.datadr, testCase.DataDir);
            testCase.verifyEqual(tt.savedr, testCase.SaveDir);
            testCase.verifyEqual(tt.filename, [dmd1; dmd2]);
            testCase.verifyEqual(double(tt.slap2_info.last_line), ...
                [100000 120000; 100000 120000]);
            testCase.verifyTrue(isfield(tt.slap2_info.ref_stack, 'Path1'));
            testCase.verifyTrue(isfield(tt.slap2_info.ref_stack, 'Path2'));
        end

        function testNoReferenceStackError(testCase)
            testCase.writeTriggeredFiles(100000, 100000);

            testCase.verifyError(@() buildTrialTableSLAP2( ...
                testCase.DataDir, testCase.SaveDir, true), ?MException);
        end

        function testUnalignedDMDsError(testCase)
            testCase.writeTriggeredFiles(100000, 150000);
            testCase.writeReferenceStacks();

            testCase.verifyError(@() buildTrialTableSLAP2( ...
                testCase.DataDir, testCase.SaveDir, true), ?MException);
        end
    end

    methods (Access = private)
        function [dmd1, dmd2] = writeTriggeredFiles(testCase, nLines1, nLines2)
            % Files are timestamped one minute apart, in trial order.
            dmd1 = arrayfun(@(k) sprintf('acq_DMD1-TRIAL%06d.dat', k), ...
                1:numel(nLines1), 'UniformOutput', false);
            dmd2 = arrayfun(@(k) sprintf('acq_DMD2-TRIAL%06d.dat', k), ...
                1:numel(nLines2), 'UniformOutput', false);
            for k = 1:numel(nLines1)
                testCase.writeDat(dmd1{k}, nLines1(k), testCase.StartTime + minutes(k));
            end
            for k = 1:numel(nLines2)
                testCase.writeDat(dmd2{k}, nLines2(k), testCase.StartTime + minutes(k));
            end
        end

        function writeDat(testCase, name, nLines, modTime)
            fn = testCase.writeText(name, sprintf('%d', nLines));
            java.io.File(fn).setLastModified(int64(posixtime(modTime)*1000));
        end

        function fn = writeText(testCase, name, contents)
            fn = fullfile(testCase.DataDir, name);
            fid = fopen(fn, 'w');
            fprintf(fid, '%s', contents);
            fclose(fid);
        end

        function writeReferenceStacks(testCase)
            refDir = fullfile(testCase.DataDir, 'reference');
            mkdir(refDir);
            for dmd = 1:2
                fid = fopen(fullfile(refDir, ...
                    sprintf('acq_DMD%d_CONFIG1-REFERENCE.tif', dmd)), 'w');
                fclose(fid);
            end
        end

        function dn = datenums(testCase, names)
            dn = zeros(1, numel(names));
            for k = 1:numel(names)
                d = dir(fullfile(testCase.DataDir, names{k}));
                dn(k) = d.datenum;
            end
        end
    end
end
