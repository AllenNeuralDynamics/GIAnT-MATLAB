function [IM, meanIM, IMc, aData, peaks, discardFrames, halfIMs] = loadAndProcessTrialAsync(mocoDr, fn, numChannels, params) %#ok<INUSD>
% Fake for SILo tests: returns the per-trial images stored in a sidecar
% <fn>.fake.mat next to the registered file instead of loading the movie
% and localizing sources. Alignment data is read from the real
% *_ALIGNMENTDATA.h5 file, as in the real function. For patchwise alignment
% the two half-trial images are the mean image plus independent noise.
fake = load(fullfile(mocoDr, [fn '.fake.mat']));
IM = [];
meanIM = fake.meanIM;
IMc = fake.actIM;
peaks = fake.peaks;
discardFrames = fake.discardFrames;
halfIMs = [];
if isfield(params, 'interTrialAlignment') && strcmpi(params.interTrialAlignment, 'patchwise')
    halfIMs = meanIM + 0.5*sqrt(max(meanIM, 0)).*randn([size(meanIM, [1 2 3]) 2]);
end

fnStemEnd = strfind(fn, '_REGISTERED') - 1;
aData = loadStructFromH5(fullfile(mocoDr, [fn(1:fnStemEnd) '_ALIGNMENTDATA.h5']));
aData.dsFac = 1;
end
