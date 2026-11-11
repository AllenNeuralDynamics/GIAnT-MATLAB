function E = processAllTrials_Async(dr, fns, fls, els, selPix, sources, discardFrames, alignData, meanAligned, motOutput, roiData, validTrials, params) %#ok<INUSD>
% Fake for SILo tests: records its inputs in dr/processAllTrials_Async_args.mat
% and returns minimal per-trial extraction results for the valid trials,
% with the fields saveExperimentSummaryH5/savePerTrialSummaryH5 read.
save(fullfile(dr, 'processAllTrials_Async_args.mat'), 'dr', 'fns', 'fls', 'els', ...
    'selPix', 'sources', 'discardFrames', 'motOutput', 'roiData', 'validTrials', 'params');

nSources = numel(sources.R);
nSelPix = nnz(any(selPix, 3));
E = cell(numel(fls), 1);
for trialIx = validTrials(:)'
    n = numel(discardFrames{trialIx});
    E{trialIx} = struct('discardFrames', discardFrames{trialIx}, ...
        'frameLines', 1:n, ...
        'global', struct('F', ones(1, n)), ...
        'footprints', ones(nSelPix, nSources), ...
        'dF', struct('ls', zeros(nSources, n, 1)), ...
        'SNR', ones(nSources, 1));
end
end
