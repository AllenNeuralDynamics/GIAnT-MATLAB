function [bytesPerTrial, channelBytes] = localizationTrialBytes(mocoDr, fns, numChannels)
%LOCALIZATIONTRIALBYTES Memory a worker needs to localize the largest of a path's trials
%   [bytesPerTrial, channelBytes] = localizationTrialBytes(mocoDr, fns, numChannels)
%   fns: the registered (downsampled) movies that loadAndProcessTrialAsync
%   reads, one per trial to localize (a cell array of names in mocoDr).
%   channelBytes: bytes of one channel of each movie as single, one per file
%   (the size of the data, not of the file: h5 movies may be compressed).
%   bytesPerTrial: what a worker needs to localize the largest movie.
%
%   Every movie is sized because trials differ in length: the first and last
%   trials of a recording are often partial, and a partial first trial can be
%   a small fraction of the others.
%
%   A worker peaks at ~4.3x one channel of its movie plus ~0.8 GB, while
%   localizeSources_vIM smooths the movie and makes its baseline. It then
%   holds four movie-sized single arrays of the activity channel (the movie,
%   its variance factors, the smoothed movie and the baseline) and its NaN
%   mask. Loading the movie's two channels peaks at ~4x. Measured on 780125
%   trials of ~2250 frames (2.3-3.0 GiB per channel, TIF and h5, Windows).
%   Trials of up to ~500 frames peak at ~5.7x: the variance of the first
%   500 valid frames copies the whole smoothed movie. Allow 6x.

channelBytes = nan(1, numel(fns));
for ix = 1:numel(fns)
    fn = fullfile(mocoDr, fns{ix});
    try
        if endsWith(fn, '.h5')
            desc = h5info(fn);
            channelBytes(ix) = 4*prod(desc.Datasets(1).Dataspace.Size)/numChannels;
        else
            dd = dir(fn);
            channelBytes(ix) = dd.bytes/numChannels;
        end
    catch
        error('localizationTrialBytes:unreadable', ...
            'Error loading registered movie: %s\nAre paths in your trial table valid?', fns{ix});
    end
end
bytesPerTrial = 6*max(channelBytes) + 2^30;
end
