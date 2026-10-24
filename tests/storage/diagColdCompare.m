function diagColdCompare(dataDir, firstRank)
%DIAGCOLDCOMPARE MATLAB vs Python concurrent reads on never-read files (one file each).
% Uses tiffs ranked firstRank.. by size (default 10) so earlier tests have not read them.
if nargin < 2, firstRank = 10; end
pool = backgroundPool; fetchOutputs(parfeval(pool, @() 1, 1));
files = dir(fullfile(dataDir, '**', '*.tif'));
files = files([files.bytes] <= 3e9 & [files.bytes] >= 5e8);
[~, o] = sort([files.bytes], 'descend');
files = files(o(firstRank:min(firstRank + 3, end)));
py = [tempname '.py'];
fid = fopen(py, 'w');
fprintf(fid, '%s\n', 'import sys, os, time', 'from concurrent.futures import ThreadPoolExecutor', ...
    'p, n = sys.argv[1], int(sys.argv[2]); size = os.path.getsize(p); e = [round(i*size/n) for i in range(n+1)]', ...
    'def part(k):', '    with open(p, "rb", buffering=0) as f:', '        f.seek(e[k]); r = e[k+1]-e[k]; got = 0', ...
    '        while r > 0:', '            c = f.read(min(r, 64000000)); got += len(c); r -= len(c)', '            if not c: break', ...
    '        return got', 't = time.time()', 'with ThreadPoolExecutor(n) as ex: tot = sum(ex.map(part, range(n)))', ...
    'print("%.1f" % (time.time() - t))');
fclose(fid);
labels = {'MATLAB dataPlanes, 8 parallel', 'Python, 8 threads', 'MATLAB single read', 'Python, 1 thread'};
for k = 1:numel(files)
    f = fullfile(files(k).folder, files(k).name); GB = files(k).bytes / 1e9;
    switch k
        case 1, r = ScanImageTiffReader(f); n = r.shape().dims(end); t = tic; d = r.dataPlanes(1:n, 8); sec = toc(t); r.close(); clear d
        case 2, [~, out] = system(sprintf('python3 "%s" "%s" 8', py, f)); sec = str2double(strtrim(out));
        case 3, r = ScanImageTiffReader(f); t = tic; d = r.data(); sec = toc(t); r.close(); clear d
        case 4, [~, out] = system(sprintf('python3 "%s" "%s" 1', py, f)); sec = str2double(strtrim(out));
    end
    fprintf('%-30s %6.1f s %6.0f MB/s  (%.2f GB, never read before)  %s\n', labels{k}, sec, GB * 1e3 / sec, GB, files(k).name);
end
delete(py);
end
