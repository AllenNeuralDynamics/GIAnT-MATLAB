function E = processAllTrials_Async(dr, fns, fls, els, selPix, sources, discardFrames, alignData, meanAligned, warp, roiData, validTrials, params) %#ok<INUSD>
% Fake for SILo tests: records its inputs in dr/processAllTrials_Async_args.mat
% (with motOutput, the offset part of each trial's map) and returns minimal
% per-trial extraction results for the valid trials, with the fields
% saveExperimentSummaryH5/savePerTrialSummaryH5 read.
motOutput = reshape(warp.affine(:,3,:), 2, []);
save(fullfile(dr, 'processAllTrials_Async_args.mat'), 'dr', 'fns', 'fls', 'els', ...
    'selPix', 'sources', 'discardFrames', 'warp', 'motOutput', 'roiData', 'validTrials', 'params');

nSources = numel(sources.R);
nSelPix = nnz(any(selPix, 3));
E = cell(numel(fls), 1);
for trialIx = validTrials(:)'
    n = numel(discardFrames{trialIx});
    E{trialIx} = struct('discardFrames', discardFrames{trialIx}, ...
        'frameLines', 1:n, ...
        'global', struct('F', ones(1, n)), ...
        'footprints', ones(nSelPix, nSources), ...
        'dF', struct('ls', fakeTraces(nSources, n, trialIx, params)), ...
        'SNR', ones(nSources, 1));
end
end

function ls = fakeTraces(nSources, n, trialIx, params)
% Zeros unless the frame rate resolves tau_s (at least 2 frames, as
% sourceKinetics needs); then transients every 40 frames in noise, decaying
% with tau_s for even-numbered sources and lasting one frame for odd ones.
ls = zeros(nSources, n, 1);
tauF = params.tau_s * params.analyzeHz;
if tauF < 2
    return
end
rs = RandStream('twister', 'Seed', trialIx);
k = exp(-(0:ceil(6 * tauF)) / tauF);
ev = zeros(1, n);
ev(20:40:n-20) = 10;
slow = conv(ev, k);
for s = 1:nSources
    if mod(s, 2) == 0
        ls(s, :) = slow(1:n);
    else
        ls(s, :) = ev;
    end
    ls(s, :) = ls(s, :) + 0.5 * randn(rs, 1, n);
end
end
