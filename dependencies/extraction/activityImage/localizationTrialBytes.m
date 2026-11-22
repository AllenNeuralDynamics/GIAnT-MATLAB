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
%   A worker peaks at ~8.3x one channel of its movie plus ~1 GB, while
%   localizeSources_vIM computes the s.d. (stdIM). It then holds eight
%   movie-sized single arrays of the activity channel and its NaN mask: the
%   movie and the variance factors that loadAndProcessTrialAsync keeps, the
%   masked copies localizeSources_vIM makes of both, the smoothed movie, the
%   baseline, the s.d. and one temporary. The other channels are freed before
%   localization (loading two channels peaks at ~4x). Measured: 26.0 GB for a
%   780125 trial of 3.04 GiB per channel (Windows); 5.7-6.7x or more on Linux
%   exports (lower bounds from total memory with 5-6 workers). Allow 9x.

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
bytesPerTrial = 9*max(channelBytes) + 2^30;
end
