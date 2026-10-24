function diagReadStackPath(dataDir)
%DIAGREADSTACKPATH Compare readStack with a direct parallel dataPlanes call on the same files.

pool = backgroundPool; f0 = parfeval(pool, @() 1, 1); fetchOutputs(f0); %start the pool before timing
files = dir(fullfile(dataDir, '**', '*.tif'));
files = files([files.bytes] <= 3e9 & [files.bytes] >= 5e8);
[~, o] = sort([files.bytes], 'descend');
files = files(o(1:min(3, end)));
for k = 1:numel(files)
    f = fullfile(files(k).folder, files(k).name); GB = files(k).bytes / 1e9;
    r = ScanImageTiffReader(f); n = r.shape().dims(end);
    t = tic; d = r.dataPlanes(1:n, 8); tDirect = toc(t); r.close(); clear d
    t = tic; [rr, sel, nPar, info, c] = openStackReader(f); tOpen = toc(t);
    t = tic; d = rr.dataPlanes(sel, nPar); tRead = toc(t); clear d c rr
    t = tic; d = readStack(f); tStack = toc(t); clear d
    fprintf('%s (%.2f GB): dataPlanes(8) %.1f s (%.0f MB/s) | openStackReader %.1f s [%s, nPar %d] + read %.1f s | readStack %.1f s (%.0f MB/s)\n', ...
        files(k).name, GB, tDirect, GB * 1e3 / tDirect, tOpen, info.strategy, nPar, tRead, tStack, GB * 1e3 / tStack);
end
end
