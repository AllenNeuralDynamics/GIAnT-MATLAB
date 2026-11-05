classdef TestXcorr2Nans < matlab.unittest.TestCase
    %TESTXCORR2NANS Unit tests for xcorr2_nans.
    %   Frames are built from a smooth random texture so that
    %   frame(j) = template(j - d) inside an ROI and NaN elsewhere; the
    %   expected motion is then d.
    %   Requires Image Processing Toolbox (strel, imdilate).
    %   Run with: runtests('tests/unit/dependencies/alignment/TestXcorr2Nans.m')

    properties
        Template
        Roi = 16:45
    end

    methods (TestClassSetup)
        function addPaths(testCase)
            import matlab.unittest.fixtures.PathFixture
            testsDir = fileparts(mfilename('fullpath'));
            repoRoot = fileparts(fileparts(fileparts(fileparts(testsDir))));
            testCase.applyFixture(PathFixture(fullfile(repoRoot, 'dependencies', 'alignment')));
        end

        function requireImageProcessingToolbox(testCase)
            testCase.assumeNotEmpty(ver('images'), ...
                'xcorr2_nans requires Image Processing Toolbox.');
        end

        function makeTemplate(testCase)
            testCase.Template = TestXcorr2Nans.makeTexture(0);
        end
    end

    methods (Test)
        function testRecoverIntegerShift(testCase)
            d = [1 -2];
            frame = testCase.shiftedFrame(testCase.Template, d);

            [motion, R] = xcorr2_nans(frame, testCase.Template, [0; 0], 3);

            testCase.verifySize(motion, [1 2]);
            testCase.verifyEqual(motion, d, 'AbsTol', 0.25);
            testCase.verifyEqual(R, 1, 'AbsTol', 1e-10);
        end

        function testRecoverSubpixelShift(testCase)
            % Averaging adjacent rows approximates a half-pixel row shift.
            frame = nan(size(testCase.Template));
            roi = testCase.Roi;
            frame(roi, roi) = (testCase.Template(roi, roi) + testCase.Template(roi - 1, roi))/2;

            motion = xcorr2_nans(frame, testCase.Template, [0; 0], 3);

            testCase.verifyEqual(motion, [0.5 0], 'AbsTol', 0.2);
        end

        function testSearchAroundShiftsCenter(testCase)
            % The true shift is outside +/-dShift of zero but within
            % +/-dShift of shiftsCenter.
            d = [5 -4];
            frame = testCase.shiftedFrame(testCase.Template, d);

            motion = xcorr2_nans(frame, testCase.Template, [4; -4], 2);

            testCase.verifyEqual(motion, d, 'AbsTol', 0.25);
        end

        function testShiftBeyondSearchRangeClampsToEdge(testCase)
            % The best match is on the edge of the search range, so no
            % subpixel refinement is applied and motion is an integer.
            frame = testCase.shiftedFrame(testCase.Template, [3 -3]);

            motion = xcorr2_nans(frame, testCase.Template, [0; 0], 2);

            testCase.verifyEqual(motion, [2 -2]);
        end

        function testInvariantToIntensityScale(testCase)
            d = [1 -2];
            frame = testCase.shiftedFrame(testCase.Template, d);

            [motion1, R1] = xcorr2_nans(frame, testCase.Template, [0; 0], 3);
            [motion2, R2] = xcorr2_nans(3*frame, testCase.Template, [0; 0], 3);

            testCase.verifyEqual(motion2, motion1, 'AbsTol', 1e-9);
            testCase.verifyEqual(R2, R1, 'AbsTol', 1e-9);
        end

        function testIgnoresFrameNaNs(testCase)
            d = [1 -2];
            frame = testCase.shiftedFrame(testCase.Template, d);
            frame(20:27, 30:36) = nan;

            [motion, R] = xcorr2_nans(frame, testCase.Template, [0; 0], 3);

            testCase.verifyEqual(motion, d, 'AbsTol', 0.25);
            testCase.verifyEqual(R, 1, 'AbsTol', 1e-10);
        end

        function testIgnoresTemplateNaNs(testCase)
            d = [1 -2];
            frame = testCase.shiftedFrame(testCase.Template, d);
            template = testCase.Template;
            template(25:30, 25:30) = nan;

            [motion, R] = xcorr2_nans(frame, template, [0; 0], 3);

            testCase.verifyFalse(any(isnan(motion)));
            testCase.verifyEqual(motion, d, 'AbsTol', 0.25);
            testCase.verifyEqual(R, 1, 'AbsTol', 1e-10);
        end

        function testRoundsDShift(testCase)
            frame = testCase.shiftedFrame(testCase.Template, [1 -2]);

            [motion1, R1] = xcorr2_nans(frame, testCase.Template, [0; 0], 3);
            [motion2, R2] = xcorr2_nans(frame, testCase.Template, [0; 0], 3.4);

            testCase.verifyEqual(motion2, motion1);
            testCase.verifyEqual(R2, R1);
        end
    end

    methods (Access = private)
        function frame = shiftedFrame(testCase, tex, d)
            % frame(j) = tex(j - d) inside the ROI, NaN elsewhere.
            roi = testCase.Roi;
            frame = nan(size(tex));
            frame(roi, roi) = tex(roi - d(1), roi - d(2));
        end
    end

    methods (Static, Access = private)
        function tex = makeTexture(seed)
            % Smooth 60x60 random texture in [10, 100].
            stream = RandStream('mt19937ar', 'Seed', seed);
            noise = randn(stream, 72);
            [x, y] = meshgrid(-6:6);
            kernel = exp(-(x.^2 + y.^2)/8);
            tex = conv2(noise, kernel/sum(kernel(:)), 'valid');
            tex = 10 + 90*(tex - min(tex(:)))/(max(tex(:)) - min(tex(:)));
        end
    end
end
