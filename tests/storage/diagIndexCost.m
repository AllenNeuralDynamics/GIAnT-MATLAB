function diagIndexCost(dataDir, firstRank)
%DIAGINDEXCOST Time opening (page index) and reading data, cold, on never-read files.
% Uses tiffs ranked firstRank.. by size (default 20) so earlier tests have not
% read them: file 1 indexed by walking the IFD chain, file 2 with the predicted
% index, file 3 the full readStack. Each file is read once (cold).
if nargin < 2, firstRank = 20; end
pool = backgroundPool; fetchOutputs(parfeval(pool, @() 1, 1));
files = dir(fullfile(dataDir, '**', '*.tif'));
files = files([files.bytes] <= 3e9 & [files.bytes] >= 5e8);
[~, o] = sort([files.bytes], 'descend');
files = files(o(firstRank:min(firstRank + 2, end)));
modes = {'walk', 'predict'};
for k = 1:min(2, numel(files))
    f = fullfile(files(k).folder, files(k).name); GB = files(k).bytes / 1e9;
    t = tic; r = ScanImageTiffReader(f, modes{k}); tOpen = toc(t); n = r.shape().dims(end);
    t = tic; d = r.dataPlanes(1:n, 8); tRead = toc(t); used = r.indexMethod(); r.close(); clear d
    fprintf('%-7s index: open %5.1f s (%s, %d pages), read data (8 parallel) %4.1f s; total %5.1f s (%.0f MB/s)  %s\n', ...
        modes{k}, tOpen, used, n, tRead, tOpen + tRead, GB * 1e3 / (tOpen + tRead), files(k).name);
end
if numel(files) >= 3
    f = fullfile(files(3).folder, files(3).name); GB = files(3).bytes / 1e9;
    t = tic; d = readStack(f); sec = toc(t); clear d
    [~, ~, ~, info, c] = openStackReader(f); clear c
    fprintf('readStack total %5.1f s (%.0f MB/s; %s, %s index)  %s\n', sec, GB * 1e3 / sec, info.strategy, info.index, files(3).name);
end
end
