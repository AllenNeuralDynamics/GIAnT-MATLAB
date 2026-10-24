function diagIndexCost(dataDir, firstRank)
%DIAGINDEXCOST Split a cold read into opening (page index) and reading data, on never-read files.
if nargin < 2, firstRank = 20; end
pool = backgroundPool; fetchOutputs(parfeval(pool, @() 1, 1));
files = dir(fullfile(dataDir, '**', '*.tif'));
files = files([files.bytes] <= 3e9 & [files.bytes] >= 5e8);
[~, o] = sort([files.bytes], 'descend');
files = files(o(firstRank:min(firstRank + 1, end)));
for k = 1:numel(files)
    f = fullfile(files(k).folder, files(k).name); GB = files(k).bytes / 1e9;
    t = tic; r = ScanImageTiffReader(f); tOpen = toc(t); n = r.shape().dims(end);
    t = tic; d = r.dataPlanes(1:n, 8); tRead = toc(t); r.close(); clear d
    fprintf('%s (%.2f GB, %d pages, never read): open/index %.1f s, read data (8 parallel) %.1f s (%.0f MB/s)\n', ...
        files(k).name, GB, n, tOpen, tRead, GB * 1e3 / tRead);
end
end
