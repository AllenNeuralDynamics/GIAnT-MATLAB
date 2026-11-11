classdef TestXcorr2NansWeighted < matlab.unittest.TestCase
    %TESTXCORR2NANSWEIGHTED Unit tests for xcorr2_nans_weighted.
    %   Frames are built from a smooth random texture so that
    %   frame(j) = template(j - d) inside an ROI and NaN elsewhere; the
    %   expected motion is then d.
    %   Run with: runtests('tests/unit/dependencies/alignment/TestXcorr2NansWeighted.m')

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

        function makeTemplate(testCase)
            testCase.Template = TestXcorr2NansWeighted.makeTexture(0);
        end
    end

    methods (Test)
        function testRecoverIntegerShift(testCase)
            d = [1 -2];
            frame = testCase.shiftedFrame(testCase.Template, d);

            [motion, R] = xcorr2_nans_weighted(frame, ones(size(frame)), ...
                testCase.Template, [0; 0], 3);

            testCase.verifySize(motion, [1 2]);
            testCase.verifyEqual(motion, d, 'AbsTol', 0.25);
            testCase.verifyEqual(R, 1, 'AbsTol', 1e-10);
        end

        function testRecoverSubpixelShift(testCase)
            % Averaging adjacent rows approximates a half-pixel row shift.
            frame = nan(size(testCase.Template));
            roi = testCase.Roi;
            frame(roi, roi) = (testCase.Template(roi, roi) + testCase.Template(roi - 1, roi))/2;

            motion = xcorr2_nans_weighted(frame, ones(size(frame)), ...
                testCase.Template, [0; 0], 3);

            testCase.verifyEqual(motion, [0.5 0], 'AbsTol', 0.2);
        end

        function testSearchAroundShiftsCenter(testCase)
            % The true shift is outside +/-dShift of zero but within
            % +/-dShift of shiftsCenter.
            d = [5 -4];
            frame = testCase.shiftedFrame(testCase.Template, d);

            motion = xcorr2_nans_weighted(frame, ones(size(frame)), ...
                testCase.Template, [4; -4], 2);

            testCase.verifyEqual(motion, d, 'AbsTol', 0.25);
        end

        function testShiftBeyondSearchRangeClampsToEdge(testCase)
            % The best match is on the edge of the search range, so no
            % subpixel refinement is applied and motion is an integer.
            frame = testCase.shiftedFrame(testCase.Template, [3 -3]);

            motion = xcorr2_nans_weighted(frame, ones(size(frame)), ...
                testCase.Template, [0; 0], 2);

            testCase.verifyEqual(motion, [2 -2]);
        end

        function testInvariantToIntensityScaleAndOffset(testCase)
            d = [1 -2];
            frame = testCase.shiftedFrame(testCase.Template, d);
            freshness = ones(size(frame));

            [motion1, R1] = xcorr2_nans_weighted(frame, freshness, testCase.Template, [0; 0], 3);
            [motion2, R2] = xcorr2_nans_weighted(2*frame + 5, freshness, testCase.Template, [0; 0], 3);

            testCase.verifyEqual(motion2, motion1, 'AbsTol', 1e-9);
            testCase.verifyEqual(R2, R1, 'AbsTol', 1e-9);
        end

        function testUniformFreshnessScalingNoEffect(testCase)
            frame = testCase.shiftedFrame(testCase.Template, [1 -2]);

            [motion1, R1] = xcorr2_nans_weighted(frame, ones(size(frame)), ...
                testCase.Template, [0; 0], 3);
            [motion2, R2] = xcorr2_nans_weighted(frame, 5*ones(size(frame)), ...
                testCase.Template, [0; 0], 3);

            testCase.verifyEqual(motion2, motion1, 'AbsTol', 1e-9);
            testCase.verifyEqual(R2, R1, 'AbsTol', 1e-9);
        end

        function testIgnoresNaNs(testCase)
            d = [1 -2];
            frame = testCase.shiftedFrame(testCase.Template, d);
            frame(20:27, 30:36) = nan;
            template = testCase.Template;
            template(25:30, 25:30) = nan;

            [motion, R] = xcorr2_nans_weighted(frame, ones(size(frame)), template, [0; 0], 3);

            testCase.verifyFalse(any(isnan(motion)));
            testCase.verifyEqual(motion, d, 'AbsTol', 0.25);
            testCase.verifyEqual(R, 1, 'AbsTol', 1e-10);
        end

        function testZeroFreshnessPixelsIgnored(testCase)
            % Replace the right half of the ROI with unrelated content and
            % give it zero weight; the match should rely on the left half.
            d = [1 -2];
            frame = testCase.shiftedFrame(testCase.Template, d);
            other = TestXcorr2NansWeighted.makeTexture(1);
            badCols = 31:45;
            frame(testCase.Roi, badCols) = other(testCase.Roi, badCols);
            freshness = ones(size(frame));
            freshness(:, badCols) = 0;

            [motion, R] = xcorr2_nans_weighted(frame, freshness, testCase.Template, [0; 0], 3);
            [~, Runweighted] = xcorr2_nans_weighted(frame, ones(size(frame)), ...
                testCase.Template, [0; 0], 3);

            testCase.verifyEqual(motion, d, 'AbsTol', 0.25);
            testCase.verifyGreaterThan(R, Runweighted);
        end

        function testRoundsDShift(testCase)
            frame = testCase.shiftedFrame(testCase.Template, [1 -2]);
            freshness = ones(size(frame));

            [motion1, R1] = xcorr2_nans_weighted(frame, freshness, testCase.Template, [0; 0], 3);
            [motion2, R2] = xcorr2_nans_weighted(frame, freshness, testCase.Template, [0; 0], 3.4);

            testCase.verifyEqual(motion2, motion1);
            testCase.verifyEqual(R2, R1);
        end

        function testMatchesPerPixelReference(testCase)
            % Identical results to the per-pixel (find/sub2ind) implementation
            % it replaced, on frames with NaN borders and holes, nonuniform
            % freshness, offset search centers and shifts past the edge.
            stream = RandStream('mt19937ar', 'Seed', 3);
            for rep = 1:12
                sz = [40 60] + randi(stream, 20, 1, 2);
                T = conv2(rand(stream, sz + 4), ones(5)/25, 'valid');
                T(rand(stream, sz) < 0.05) = nan;
                F = circshift(T, randi(stream, [-4 4], 1, 2)) + 0.05*randn(stream, sz);
                F(:, 1:randi(stream, 10)) = nan;
                F(rand(stream, sz) < 0.1) = nan;
                fresh = 1 + 3*rand(stream, sz);
                center = randi(stream, [-3 3], 2, 1);
                dShift = randi(stream, [2 6]);

                [motion, R] = xcorr2_nans_weighted(F, fresh, T, center, dShift);
                [motionRef, Rref] = perPixelReference(F, fresh, T, center, dShift);

                testCase.verifyEqual(motion, motionRef);
                testCase.verifyEqual(R, Rref);
            end
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

function [motion, R] = perPixelReference(frame, freshness, template, shiftsCenter, dShift)
% xcorr2_nans_weighted before it took overlapping blocks per shift
sz = size(template);
dShift = round(dShift(1));
template = circshift(template, shiftsCenter);
[Fr0, Fc0] = find(~isnan(frame));
shifts = -dShift:dShift;
C = nan(length(shifts), length(shifts));
for drix = 1:length(shifts)
    for dcix = 1:length(shifts)
        Tr = Fr0 + shifts(drix);
        Tc = Fc0 + shifts(dcix);
        valid = Tr>=1 & Tr<=sz(1) & Tc>=1 & Tc<=sz(2);
        Fr = Fr0(valid); Fc = Fc0(valid); Tr = Tr(valid); Tc = Tc(valid);
        T = template(sub2ind(sz, Tr, Tc));
        sel = ~isnan(T); T = T(sel);
        indsF = sub2ind(sz, Fr(sel), Fc(sel));
        F = frame(indsF);
        Ff = freshness(indsF);
        sFT = sum(Ff.*(F-mean(F)).*(T-mean(T)));
        sT = mean((T-mean(T)).^2);
        sF = sum(Ff.*(F-mean(F)).^2);
        C(drix,dcix) = sFT./sqrt(sT.*sF.*sum(Ff));
    end
end
[maxval, I] = max(C(:));
[rr, cc] = ind2sub(size(C), I);
R = maxval;
if rr>1 && rr<length(shifts) && cc>1 && cc<length(shifts)
    ratioR = min(1e6,(C(rr,cc) - C(rr-1,cc))/(C(rr,cc) - C(rr+1,cc)));
    dR = (1-ratioR)/(1+ratioR)/2;
    ratioC = min(1e6,(C(rr,cc) - C(rr,cc-1))/(C(rr,cc) - C(rr,cc+1)));
    dC = (1-ratioC)/(1+ratioC)/2;
    motion = -shiftsCenter' + [shifts(rr)-dR shifts(cc)-dC];
else
    motion = -shiftsCenter' + [shifts(rr) shifts(cc)];
end
motion = -motion;
end
