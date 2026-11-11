classdef TestLoadStructFromH5 < matlab.unittest.TestCase
    %TESTLOADSTRUCTFROMH5 Unit tests for loadStructFromH5's skipDatasets option.
    %   SILo and verifyFiles read alignment files only for their metadata,
    %   and skip /slap2/varFacDS, a variance movie as large as the
    %   registered movie.
    %   Run with: runtests('tests/unit/dependencies/TestLoadStructFromH5.m')

    properties
        Fn
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
        function writeAlignmentLikeFile(testCase)
            import matlab.unittest.fixtures.TemporaryFolderFixture
            tmp = testCase.applyFixture(TemporaryFolderFixture);
            testCase.Fn = fullfile(tmp.Folder, 'E1T1DMD1_ALIGNMENTDATA.h5');
            s.numChannels = 2;
            s.alignHz = 80;
            s.recNegErr = (1:10)';
            s.slap2.Z = 12.5; % old name for Z_depths
            s.slap2.varFacDS = rand(4, 5, 6, 'single');
            s.slap2.viewR = (1:4)';
            saveStructToH5(s, testCase.Fn);
        end
    end

    methods (Test)
        function testReadsEverythingByDefault(testCase)
            s = loadStructFromH5(testCase.Fn);
            testCase.verifyTrue(isfield(s.slap2, 'varFacDS'));
            testCase.verifySize(s.slap2.varFacDS, [4 5 6]);
        end

        function testSkipsListedDataset(testCase)
            full = loadStructFromH5(testCase.Fn);
            s = loadStructFromH5(testCase.Fn, '/slap2/varFacDS');
            testCase.verifyFalse(isfield(s.slap2, 'varFacDS'));
            % everything else is read as without skipping
            full.slap2 = rmfield(full.slap2, 'varFacDS');
            testCase.verifyEqual(s, full);
        end

        function testKeepsZDepthsNormalization(testCase)
            s = loadStructFromH5(testCase.Fn, {'/slap2/varFacDS'});
            testCase.verifyEqual(s.slap2.Z_depths, 12.5);
            testCase.verifyFalse(isfield(s.slap2, 'Z'));
        end

        function testSkipMatchesFullPathOnly(testCase)
            % a dataset of the same name in another group is still read
            s = loadStructFromH5(testCase.Fn, '/varFacDS');
            testCase.verifyTrue(isfield(s.slap2, 'varFacDS'));
        end
    end
end
