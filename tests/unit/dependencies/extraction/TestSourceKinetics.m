classdef TestSourceKinetics < matlab.unittest.TestCase
    %TESTSOURCEKINETICS Unit tests for sourceKinetics and rejectFastSources.
    %   Synthetic least-squares traces with transients that decay like the
    %   indicator (tau_s), one-frame transients, or noise only. The kinetics
    %   score should be about 1, about 0 and NaN for these, and
    %   rejectFastSources should remove only the one-frame sources from the
    %   trials' extraction results and the seeds.
    %   Run with: runtests('tests/unit/dependencies/extraction/TestSourceKinetics.m')

    properties
        Fs = 200
        Tau_s = 0.024
        NumTrials = 3
        NumFrames = 2000
        Kinds = {'indicator', 'oneFrame', 'noise'}
    end

    methods (TestClassSetup)
        function addPaths(testCase)
            import matlab.unittest.fixtures.PathFixture
            testsDir = fileparts(mfilename('fullpath'));
            repoRoot = fileparts(fileparts(fileparts(fileparts(testsDir))));
            testCase.applyFixture(PathFixture(fullfile(repoRoot, 'dependencies', 'extraction')));
        end
    end

    methods (Test)
        function testScoresByKinetics(testCase)
            E = testCase.makeTrials();
            [score, nEvents] = sourceKinetics(E, testCase.params());
            testCase.verifyEqual(score(1), 1, 'AbsTol', 0.25);
            testCase.verifyEqual(score(2), 0, 'AbsTol', 0.25);
            testCase.verifyTrue(isnan(score(3)));
            testCase.verifyGreaterThanOrEqual(nEvents(1:2), 50);
            testCase.verifyLessThan(nEvents(3), 5);
        end

        function testDiscardedFramesAndEmptyTrials(testCase)
            E = testCase.makeTrials();
            discard = false(1, testCase.NumFrames);
            discard(1:10:end) = true;
            for t = 1:numel(E)
                E{t}.dF.ls(:, discard, :) = nan;
            end
            E = [{[]}; E; {struct('discardFrames', discard)}];
            score = sourceKinetics(E, testCase.params());
            testCase.verifyEqual(score(1), 1, 'AbsTol', 0.25);
            testCase.verifyEqual(score(2), 0, 'AbsTol', 0.25);
        end

        function testUnresolvedKineticsGiveNaN(testCase)
            % at 20 Hz the 24 ms decay is under one frame
            E = testCase.makeTrials();
            p = testCase.params();
            p.analyzeHz = 20;
            [score, nEvents] = sourceKinetics(E, p);
            testCase.verifyTrue(all(isnan(score)));
            testCase.verifyEqual(nEvents, zeros(3, 1));
        end

        function testNoTracesGiveEmptyScores(testCase)
            [score, nEvents] = sourceKinetics({[], struct('SNR', 1)}, testCase.params());
            testCase.verifyEmpty(score);
            testCase.verifyEmpty(nEvents);
        end

        function testRemovesOneFrameSources(testCase)
            E = testCase.makeTrials();
            sources = struct('R', [11; 12; 13], 'C', [21; 22; 23], 'V', [1; 2; 3]);
            p = testCase.params();
            p.minKineticScore = 0.4;
            [E2, s2, kin] = rejectFastSources(E, sources, p);

            keep = [1; 3]; % the one-frame source goes; the unscored one stays
            testCase.verifyEqual(s2.R, sources.R(keep));
            testCase.verifyEqual(s2.C, sources.C(keep));
            testCase.verifyEqual(s2.V, sources.V(keep));
            for t = 1:numel(E)
                testCase.verifyEqual(E2{t}.footprints, E{t}.footprints(:, keep));
                testCase.verifyEqual(E2{t}.dF.ls, E{t}.dF.ls(keep, :, :));
                testCase.verifyEqual(E2{t}.dF.denoised, E{t}.dF.denoised(keep, :, :));
                testCase.verifyEqual(E2{t}.dF.events, E{t}.dF.events(keep, :, :));
                testCase.verifyEqual(E2{t}.F0, E{t}.F0(keep, :, :));
                testCase.verifyEqual(E2{t}.SNR, E{t}.SNR(keep));
                testCase.verifyEqual(E2{t}.discardFrames, E{t}.discardFrames);
            end
            testCase.verifySize(kin.score, [2 1]);
            testCase.verifyTrue(isnan(kin.score(2)));
            testCase.verifyEqual(kin.rejected.R, 12);
            testCase.verifyEqual(kin.rejected.C, 22);
            testCase.verifyLessThan(kin.rejected.score, 0.4);
        end

        function testZeroCutoffKeepsEverySource(testCase)
            E = testCase.makeTrials();
            sources = struct('R', [11; 12; 13], 'C', [21; 22; 23], 'V', [1; 2; 3]);
            p = testCase.params();
            p.minKineticScore = 0;
            [E2, s2, kin] = rejectFastSources(E, sources, p);
            testCase.verifyEqual(s2, sources);
            testCase.verifyEqual(E2, E);
            testCase.verifySize(kin.score, [3 1]);
            testCase.verifyEmpty(kin.rejected.R);
        end
    end

    methods
        function p = params(testCase)
            p = struct('tau_s', testCase.Tau_s, 'analyzeHz', testCase.Fs);
        end

        function E = makeTrials(testCase)
            % one source per kind, in unit Gaussian noise, two channels:
            % transients every 100 frames (jittered), amplitude 10 decaying
            % with tau_s, or amplitude 25 for one frame (the fast transients
            % seen on SLAP2 are large single-frame excursions)
            rs = RandStream('twister', 'Seed', 7);
            tauF = testCase.Tau_s * testCase.Fs;
            k = exp(-(0:ceil(10 * tauF)) / tauF);
            nS = numel(testCase.Kinds);
            T = testCase.NumFrames;
            E = cell(testCase.NumTrials, 1);
            for t = 1:testCase.NumTrials
                ls = randn(rs, nS, T, 2);
                for s = 1:nS
                    ev = zeros(1, T);
                    on = 50:100:T-50;
                    ev(on + randi(rs, [-10 10], size(on))) = 1;
                    switch testCase.Kinds{s}
                        case 'indicator'
                            tr = conv(10 * ev, k);
                            ls(s, :, 1) = ls(s, :, 1) + tr(1:T);
                        case 'oneFrame'
                            ls(s, :, 1) = ls(s, :, 1) + 25 * ev;
                    end
                end
                E{t} = struct('discardFrames', false(1, T), ...
                    'footprints', rand(rs, 50, nS), ...
                    'dF', struct('ls', ls, 'denoised', ls + 1, 'events', ls + 2), ...
                    'F0', 10 + ls, ...
                    'SNR', [5; 6; 7]);
            end
        end
    end
end
