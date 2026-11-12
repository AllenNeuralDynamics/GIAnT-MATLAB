classdef TestSaveStructToH5 < matlab.unittest.TestCase
    %TESTSAVESTRUCTTOH5 Unit tests for saveStructToH5.
    %   Every GIAnT output (trial tables, alignment data, summaries) is
    %   written by saveStructToH5 and read by loadStructFromH5 or outside
    %   MATLAB, so the file layout must stay as h5create/h5write made it.
    %   Run with: runtests('tests/unit/dependencies/TestSaveStructToH5.m')

    properties
        Fn
    end

    methods (TestClassSetup)
        function addPaths(testCase)
            import matlab.unittest.fixtures.PathFixture
            testsDir = fileparts(mfilename('fullpath'));
            repoRoot = fileparts(fileparts(fileparts(testsDir)));
            testCase.applyFixture(PathFixture(fullfile(repoRoot, 'dependencies', 'io')));
        end
    end

    methods (TestMethodSetup)
        function makeFilename(testCase)
            import matlab.unittest.fixtures.TemporaryFolderFixture
            tmp = testCase.applyFixture(TemporaryFolderFixture);
            testCase.Fn = fullfile(tmp.Folder, 'out.h5');
        end
    end

    methods (Test)
        function testRoundTripsSupportedValues(testCase)
            s.d = reshape(1:24, 2, 3, 4);
            s.f = single([1.5 NaN Inf]);
            s.i = int16([-3; 4]);
            s.u = uint64(2^60);
            s.name = 'trial 1';
            s.unit = "µm";
            s.labels = ["a" "bc"; "d" ""];
            s.files = {'x.tif', '', "y.tif"};
            s.group.nested.value = 7;

            saveStructToH5(s, testCase.Fn);
            r = loadStructFromH5(testCase.Fn);

            for f = {'d', 'f', 'i', 'u', 'name'}
                testCase.verifyEqual(r.(f{1}), s.(f{1}), f{1});
            end
            testCase.verifyEqual(r.unit, char(s.unit));
            testCase.verifyEqual(r.labels, cellstr(s.labels));
            testCase.verifyEqual(r.files, {'x.tif', '', 'y.tif'});
            testCase.verifyEqual(r.group.nested.value, 7);
        end

        function testWritesLogicalsAsInt8(testCase)
            saveStructToH5(struct('mask', [true false; false true]), testCase.Fn);

            testCase.verifyEqual(h5read(testCase.Fn, '/mask'), int8([1 0; 0 1]));
        end

        function testKeepsHighLevelWriterLayout(testCase)
            % variable-length UTF-8 strings, native types, contiguous datasets
            saveStructToH5(struct('s', "abc", 'x', zeros(3, 2, 'single')), testCase.Fn);

            str = h5info(testCase.Fn, '/s');
            testCase.verifyEqual(str.Datatype.Class, 'H5T_STRING');
            testCase.verifyEqual(str.Datatype.Type.Length, 'H5T_VARIABLE');
            testCase.verifyEqual(str.Datatype.Type.CharacterSet, 'H5T_CSET_UTF8');
            testCase.verifyEqual(str.Dataspace.Size, [1 1]);
            x = h5info(testCase.Fn, '/x');
            testCase.verifyEqual(x.Datatype.Type, 'H5T_IEEE_F32LE');
            testCase.verifyEqual(x.Dataspace.Size, [3 2]);
            testCase.verifyEmpty(x.ChunkSize);
        end

        function testWritesRowMajorFlagOnce(testCase)
            % a struct read back with loadStructFromH5 carries row_major
            saveStructToH5(struct('x', 1, 'row_major', uint8(1)), testCase.Fn);

            info = h5info(testCase.Fn);
            testCase.verifyEqual(nnz(strcmp({info.Datasets.Name}, 'row_major')), 1);
            testCase.verifyEqual(h5read(testCase.Fn, '/row_major'), uint8(0));
        end

        function testSkipsEmptyValuesAndTheirGroups(testCase)
            s.keep = 1;
            s.none = [];
            s.noText = {};
            s.emptyGroup = struct();
            s.groupOfEmpties.a = [];

            saveStructToH5(s, testCase.Fn);

            info = h5info(testCase.Fn);
            testCase.verifyEmpty(info.Groups);
            testCase.verifyEqual(sort({info.Datasets.Name}), {'keep', 'row_major'});
        end

        function testOverwritesExistingFile(testCase)
            saveStructToH5(struct('old', 1), testCase.Fn);

            saveStructToH5(struct('new', 2), testCase.Fn);

            r = loadStructFromH5(testCase.Fn);
            testCase.verifyFalse(isfield(r, 'old'));
            testCase.verifyEqual(r.new, 2);
        end

        function testRejectsUnsupportedValuesBeforeWriting(testCase)
            saveStructToH5(struct('previous', 1), testCase.Fn);
            bad = {
                sparse([1 0; 0 2]),    'saveStructToH5:UnsupportedNumeric'
                [1+2i 3],              'saveStructToH5:UnsupportedNumeric'
                ['ab'; 'cd'],          'saveStructToH5:UnsupportedChar'
                {{1, 2}},              'saveStructToH5:UnsupportedCell'
                @sin,                  'saveStructToH5:UnsupportedType'
                struct('a', {1, 2}),   'saveStructToH5:NonScalarStruct'
                };
            for k = 1:size(bad, 1)
                s = struct('ok', 1);
                s.bad = bad{k, 1};
                if iscell(s.bad), s.bad = s.bad{1}; end

                testCase.verifyError(@() saveStructToH5(s, testCase.Fn), bad{k, 2});

                r = loadStructFromH5(testCase.Fn);
                testCase.verifyEqual(r.previous, 1, ...
                    sprintf('case %d: the existing file must be left as it was', k));
            end
        end
    end
end
