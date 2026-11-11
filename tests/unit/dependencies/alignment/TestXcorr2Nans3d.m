classdef TestXcorr2Nans3d < matlab.unittest.TestCase
    %TESTXCORR2NANS3D Unit tests for xcorr2_nans3d.
    %   The template stack is built from smooth random textures; the frame
    %   is one plane (X) shifted so that frame(j) = X(j - d) inside an ROI
    %   and NaN elsewhere. The expected motion is then [d planeIndexOfX].
    %   Requires Image Processing Toolbox (strel, imdilate).
    %   Run with: runtests('tests/unit/dependencies/alignment/TestXcorr2Nans3d.m')

    properties
        X   % plane the frame is taken from
        A
        B
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
                'xcorr2_nans3d requires Image Processing Toolbox.');
        end

        function makePlanes(testCase)
            testCase.X = TestXcorr2Nans3d.makeTexture(0);
            testCase.A = TestXcorr2Nans3d.makeTexture(1);
            testCase.B = TestXcorr2Nans3d.makeTexture(2);
        end
    end

    methods (Test)
        function testFindsPlaneAndShift(testCase)
            % Symmetric neighbours around X make the z refinement exactly 0.
            d = [1 -2];
            stack = cat(3, testCase.A, testCase.B, testCase.X, testCase.B, testCase.A);
            frame = testCase.shiftedFrame(testCase.X, d);

            [motion, R] = xcorr2_nans3d(frame, stack, [0; 0], 3);

            testCase.verifySize(motion, [1 3]);
            testCase.verifyEqual(motion(1:2), d, 'AbsTol', 0.25);
            testCase.verifyEqual(motion(3), 3);
            testCase.verifyEqual(R, 1, 'AbsTol', 1e-10);
        end

        function testEdgePlaneNoSubpixel(testCase)
            % The best plane is the first one, so no subpixel refinement is
            % applied on any axis.
            stack = cat(3, testCase.X, testCase.B, testCase.A);
            frame = testCase.shiftedFrame(testCase.X, [1 -2]);

            motion = xcorr2_nans3d(frame, stack, [0; 0], 3);

            testCase.verifyEqual(motion, [1 -2 1]);
        end

        function testSearchAroundShiftsCenter(testCase)
            % The true shift is outside +/-dShift of zero but within
            % +/-dShift of shiftsCenter.
            d = [5 -4];
            stack = cat(3, testCase.A, testCase.B, testCase.X, testCase.B, testCase.A);
            frame = testCase.shiftedFrame(testCase.X, d);

            motion = xcorr2_nans3d(frame, stack, [4; -4], 2);

            testCase.verifyEqual(motion, [d 3], 'AbsTol', 0.25);
        end

        function testIgnoresTemplateNaNsInAnyPlane(testCase)
            % A NaN block in one plane is excluded from every plane.
            d = [1 -2];
            planeA = testCase.A;
            planeA(25:30, 25:30) = nan;
            stack = cat(3, planeA, testCase.B, testCase.X, testCase.B, testCase.A);
            frame = testCase.shiftedFrame(testCase.X, d);

            [motion, R] = xcorr2_nans3d(frame, stack, [0; 0], 3);

            testCase.verifyFalse(any(isnan(motion)));
            testCase.verifyEqual(motion(1:2), d, 'AbsTol', 0.25);
            testCase.verifyEqual(motion(3), 3);
            testCase.verifyEqual(R, 1, 'AbsTol', 1e-10);
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
