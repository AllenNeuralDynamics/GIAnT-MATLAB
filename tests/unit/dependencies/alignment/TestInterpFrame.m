classdef TestInterpFrame < matlab.unittest.TestCase
    %TESTINTERPFRAME Unit tests for interpFrame.
    %   interpFrame bilinearly samples M1 on the grid viewR x viewC, with
    %   each of the four neighbours weighted by its freshness, and returns
    %   a variance multiplier vIM = 1/sum(weight .* freshness).
    %   As in MultiRoiRegistration, viewC is a row and viewR a column.
    %   Run with: runtests('tests/unit/dependencies/alignment/TestInterpFrame.m')

    methods (TestClassSetup)
        function addPaths(testCase)
            import matlab.unittest.fixtures.PathFixture
            testsDir = fileparts(mfilename('fullpath'));
            repoRoot = fileparts(fileparts(fileparts(fileparts(testsDir))));
            testCase.applyFixture(PathFixture(fullfile(repoRoot, 'dependencies', 'alignment')));
        end
    end

    methods (Test)
        function testIntegerCoordinatesReturnPixels(testCase)
            M1 = magic(6);
            freshness = reshape(1:36, 6, 6);
            viewR = (2:4)';
            viewC = 3:6;

            [IM, vIM] = interpFrame(M1, viewC, viewR, freshness);

            testCase.verifyEqual(IM, M1(viewR, viewC), 'AbsTol', 1e-12);
            testCase.verifyEqual(vIM, 1./freshness(viewR, viewC), 'AbsTol', 1e-12);
        end

        function testBilinearInterpolation(testCase)
            % Bilinear interpolation reproduces a bilinear surface exactly.
            [C, R] = meshgrid(1:6);
            M1 = R + 10*C + 0.5*R.*C;
            viewR = (2:4)' + 0.3;
            viewC = (1:4) + 0.6;

            [IM, vIM] = interpFrame(M1, viewC, viewR, ones(6));

            [Cq, Rq] = meshgrid(viewC, viewR);
            testCase.verifySize(IM, [3 4]);
            testCase.verifyEqual(IM, Rq + 10*Cq + 0.5*Rq.*Cq, 'AbsTol', 1e-10);
            testCase.verifyEqual(vIM, ones(3, 4), 'AbsTol', 1e-12);
        end

        function testFreshnessWeightsNeighbours(testCase)
            % Halfway between 0 (freshness 1) and 10 (freshness 3):
            % IM = (0.5*1*0 + 0.5*3*10)/(0.5*1 + 0.5*3), vIM = 1/2.
            [IM, vIM] = interpFrame([0 10], 1.5, 1, [1 3]);

            testCase.verifyEqual(IM, 7.5, 'AbsTol', 1e-12);
            testCase.verifyEqual(vIM, 0.5, 'AbsTol', 1e-12);
        end

        function testZeroFreshnessNeighbourIgnored(testCase)
            [IM, vIM] = interpFrame([0 10], 1.5, 1, [0 1]);

            testCase.verifyEqual(IM, 10, 'AbsTol', 1e-12);
            testCase.verifyEqual(vIM, 2, 'AbsTol', 1e-12);
        end

        function testZeroFreshnessPixelIsNaN(testCase)
            M1 = magic(3);
            freshness = ones(3);
            freshness(2, 2) = 0;

            [IM, vIM] = interpFrame(M1, 1:3, (1:3)', freshness);

            expectedIM = M1;
            expectedIM(2, 2) = nan;
            expectedV = ones(3);
            expectedV(2, 2) = nan; % 1/0 is replaced with NaN
            testCase.verifyEqual(IM, expectedIM, 'AbsTol', 1e-12);
            testCase.verifyEqual(vIM, expectedV, 'AbsTol', 1e-12);
        end

        function testFreshnessScaling(testCase)
            % Scaling freshness uniformly leaves IM unchanged and scales vIM
            % by the inverse.
            M1 = magic(6);
            viewR = (2:4)' + 0.3;
            viewC = (1:4) + 0.6;

            [IM1, vIM1] = interpFrame(M1, viewC, viewR, ones(6));
            [IM4, vIM4] = interpFrame(M1, viewC, viewR, 4*ones(6));

            testCase.verifyEqual(IM4, IM1, 'AbsTol', 1e-10);
            testCase.verifyEqual(vIM4, vIM1/4, 'AbsTol', 1e-12);
        end

        function testOutOfBoundsIsNaN(testCase)
            M1 = magic(4);
            viewR = (0:2)';
            viewC = 3:5;

            [IM, vIM] = interpFrame(M1, viewC, viewR, ones(4));

            % Row 0 and column 5 are outside the 4x4 image.
            expectedIM = nan(3);
            expectedIM(2:3, 1:2) = M1(1:2, 3:4);
            expectedV = nan(3);
            expectedV(2:3, 1:2) = 1;
            testCase.verifyEqual(IM, expectedIM, 'AbsTol', 1e-12);
            testCase.verifyEqual(vIM, expectedV, 'AbsTol', 1e-12);
        end

        function testFractionalCoordinateNeedsBothNeighbours(testCase)
            % Column 4.5 needs column 5, which is outside the 4x4 image.
            M1 = magic(4);

            IM = interpFrame(M1, [3.5 4.5], (1:2)', ones(4));

            testCase.verifyEqual(IM(:, 1), (M1(1:2, 3) + M1(1:2, 4))/2, 'AbsTol', 1e-12);
            testCase.verifyTrue(all(isnan(IM(:, 2))));
        end

        function testNaNPixelPropagatesToNeighbours(testCase)
            M1 = magic(4);
            M1(2, 2) = nan;

            IM = interpFrame(M1, (1:3) + 0.5, (1:3)', ones(4));

            % Columns 1.5 and 2.5 both use column 2.
            expected = (M1(1:3, 1:3) + M1(1:3, 2:4))/2;
            testCase.verifyEqual(IM, expected, 'AbsTol', 1e-12);
            testCase.verifyTrue(all(isnan(IM(2, 1:2))));
        end

        function testUsesFractionOfFirstCoordinate(testCase)
            % The weights come from viewC(1) and viewR(1) only, so callers
            % must pass grids with a constant fractional offset. Here 2.75
            % is treated as 2.25.
            IM = interpFrame([0 10 20 30], [1.25 2.75], 1, ones(1, 4));

            testCase.verifyEqual(IM, [2.5 12.5], 'AbsTol', 1e-12);
        end
    end
end
