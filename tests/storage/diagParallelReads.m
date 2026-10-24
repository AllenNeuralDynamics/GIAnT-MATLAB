function diagParallelReads(dataDir)
%DIAGPARALLELREADS How well concurrent reads scale on dataDir's storage, from MATLAB.
% Each configuration reads a different (cold) tiff <= 3 GB under dataDir.

p = backgroundPool;
fprintf('backgroundPool workers %d, Parallel Computing Toolbox %d, cores %d, storage %s\n', ...
    p.NumWorkers, ~isempty(ver('parallel')), feature('numcores'), storageClass(fullfile(dataDir, 'x.tif')));
files = dir(fullfile(dataDir, '**', '*.tif'));
files = files([files.bytes] <= 3e9 & [files.bytes] >= 5e8);
[~, o] = sort([files.bytes], 'descend');
files = files(o);
configs = [1 2 4 8 16];
for k = 1:min(numel(configs), numel(files))
    f = fullfile(files(k).folder, files(k).name);
    r = ScanImageTiffReader(f); n = r.shape().dims(end);
    lastwarn('');
    t = tic; d = r.dataPlanes(1:n, configs(k)); sec = toc(t); r.close(); clear d
    [~, id] = lastwarn;
    fprintf('dataPlanes, %2d concurrent reads: %6.1f s %6.0f MB/s (%.2f GB) %s\n', configs(k), sec, ...
        files(k).bytes / 1e6 / sec, files(k).bytes / 1e9, id);
end
k = numel(configs) + 1;
if numel(files) >= k && ~ispc
    f = fullfile(files(k).folder, files(k).name);
    sz = files(k).bytes; nProc = 8; blk = 4 * 2^20; nBlk = ceil(sz / blk); per = ceil(nBlk / nProc);
    tmp = [tempname(getenv('GIANT_SCRATCH')) '.tif'];
    cmd = '';
    for j = 0:nProc - 1
        cmd = [cmd sprintf('dd if="%s" of="%s" bs=%d skip=%d seek=%d count=%d conv=notrunc status=none & ', ...
            f, tmp, blk, j * per, j * per, per)]; %#ok<AGROW>
    end
    t = tic; system([cmd 'wait']); tCopy = toc(t);
    r = ScanImageTiffReader(tmp); d = r.data(); r.close(); clear d; sec = toc(t); delete(tmp);
    fprintf('8 dd processes to scratch + read: %6.1f s %6.0f MB/s (copy %.1f s; %.2f GB)\n', sec, sz / 1e6 / sec, tCopy, sz / 1e9);
end
end
