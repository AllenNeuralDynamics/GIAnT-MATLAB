classdef TestSetParams < matlab.unittest.TestCase
    %TESTSETPARAMS Unit tests for setParams' modality-dependent defaults.
    %   SILo's peak thresholds (peakth, wmfPeakth) and kinetics cutoff
    %   (minKineticScore) default to empty, and setParams fills them for the
    %   recording's modality: 3, 3 and 0.4 for SLAP2, 8, 5 and 0 otherwise.
    %   Values the caller sets are kept.
    %   Run with: runtests('tests/unit/dependencies/TestSetParams.m')

    methods (TestClassSetup)
        function addPaths(testCase)
            import matlab.unittest.fixtures.PathFixture
            testsDir = fileparts(mfilename('fullpath'));
            repoRoot = fileparts(fileparts(fileparts(testsDir)));
            testCase.applyFixture(PathFixture(fullfile(repoRoot, 'dependencies', 'gui')));
        end
    end

    methods (Test)
        function testNonSlap2Defaults(testCase)
            p = setParams('SILo', struct('isSLAP2', false));
            testCase.verifyEqual(p.peakth, 8);
            testCase.verifyEqual(p.wmfPeakth, 5);
            testCase.verifyEqual(p.minKineticScore, 0);
        end

        function testSlap2Defaults(testCase)
            p = setParams('SILo', struct('isSLAP2', true));
            testCase.verifyEqual(p.peakth, 3);
            testCase.verifyEqual(p.wmfPeakth, 3);
            testCase.verifyEqual(p.minKineticScore, 0.4);
        end

        function testLegacyMicroscopeField(testCase)
            p = setParams('SILo', struct('microscope', 'SLAP2'));
            testCase.verifyTrue(p.isSLAP2);
            testCase.verifyEqual(p.peakth, 3);
            testCase.verifyEqual(p.wmfPeakth, 3);
        end

        function testCallerValuesKept(testCase)
            p = setParams('SILo', struct('isSLAP2', true, 'peakth', 8, 'wmfPeakth', 5));
            testCase.verifyEqual(p.peakth, 8);
            testCase.verifyEqual(p.wmfPeakth, 5);
            p = setParams('SILo', struct('isSLAP2', false, 'peakth', 3, 'wmfPeakth', 4));
            testCase.verifyEqual(p.peakth, 3);
            testCase.verifyEqual(p.wmfPeakth, 4);
        end

        function testJsonNullGetsDefault(testCase)
            % SILo accepts JSON parameters; null decodes to []
            p = setParams('SILo', jsondecode('{"isSLAP2":true,"peakth":null,"wmfPeakth":4}'));
            testCase.verifyEqual(p.peakth, 3);
            testCase.verifyEqual(p.wmfPeakth, 4);
        end

        function testDefaultsWithoutDisplay(testCase)
            % with no arguments and nobody to answer the GUI, the non-SLAP2 defaults
            testCase.assumeFalse(isInteractiveSession(), 'Would open the parameter GUI.');
            p = testCase.verifyWarning(@() setParams('SILo'), 'setParams:NonInteractive');
            testCase.verifyFalse(p.isSLAP2);
            testCase.verifyEqual(p.peakth, 8);
            testCase.verifyEqual(p.wmfPeakth, 5);
        end

        function testOtherFunctionsUnchanged(testCase)
            p = setParams('MultiRoiRegistration', struct('nWorkers', 1));
            testCase.verifyFalse(isfield(p, 'peakth'));
            testCase.verifyFalse(isfield(p, 'wmfPeakth'));
            testCase.verifyFalse(isfield(p, 'minKineticScore'));
        end

        function testKineticCutoffKept(testCase)
            p = setParams('SILo', struct('isSLAP2', true, 'minKineticScore', 0));
            testCase.verifyEqual(p.minKineticScore, 0);
            p = setParams('SILo', jsondecode('{"isSLAP2":true,"minKineticScore":null}'));
            testCase.verifyEqual(p.minKineticScore, 0.4);
        end
    end
end
