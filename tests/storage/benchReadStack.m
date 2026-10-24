function results = benchReadStack(dataDir, maxGB)
%BENCHREADSTACK Compare the previous tiff read with readStack, cold, on dataDir's storage.
%
% results = benchReadStack(dataDir, maxGB)
%
% Uses the 4 largest tiffs <= maxGB (default 3) under dataDir, one read per
% file so every read is cold (provided none of them has been read on this
% machine recently): previous, readStack, readStack, previous, so neither
% method always gets the larger files. "previous" is the single direct read
% GIAnT used before readStack (ScanImageTiffReader.data). Also times reading
% 1 page in 10 with readStack.

if nargin < 2 || isempty(maxGB)
    maxGB = 3;
end
pool = backgroundPool; fetchOutputs(parfeval(pool, @() 1, 1)); %start the thread pool before timing (a one-time cost per MATLAB session)
cls = storageClass(fullfile(dataDir, 'x.tif'));
fprintf('Storage of %s: %s\n', dataDir, cls);
files = dir(fullfile(dataDir, '**', '*.tif'));
files = files([files.bytes] <= maxGB * 1e9);
[~, o] = sort([files.bytes], 'descend');
files = files(o);
if numel(files) < 5
    error('benchReadStack:FewFiles', 'Need 5 tiffs <= %g GB under %s; found %d.', maxGB, dataDir, numel(files));
end

methods = {'previous', 'readStack', 'readStack', 'previous'};
results = struct('method', {}, 'file', {}, 'GB', {}, 'seconds', {}, 'MBps', {}, 'strategy', {});
for k = 1:4
    fn = fullfile(files(k).folder, files(k).name);
    t = tic;
    if strcmp(methods{k}, 'previous')
        r = ScanImageTiffReader(fn); data = r.data(); r.close(); strategy = 'single read';
        sec = toc(t);
    else
        data = readStack(fn); %data only, like the previous read (descriptions are timed below)
        sec = toc(t);
        [~, ~, ~, info, c] = openStackReader(fn); strategy = info.strategy; clear c
    end
    results(end+1) = struct('method', methods{k}, 'file', files(k).name, 'GB', files(k).bytes / 1e9, ...
        'seconds', sec, 'MBps', files(k).bytes / 1e6 / sec, 'strategy', strategy); %#ok<AGROW>
    fprintf('%-10s %6.1f s %7.0f MB/s  (%.2f GB, %s)  %s\n', methods{k}, sec, files(k).bytes / 1e6 / sec, ...
        files(k).bytes / 1e9, strategy, files(k).name);
    clear data
end

fn = fullfile(files(5).folder, files(5).name);
r = ScanImageTiffReader(fn); n = r.shape().dims(end);
t = tic; r.descriptions(); tDesc = toc(t); r.close();
fprintf('page descriptions (all %d pages, as ScanImageTiffWrapper callers request them): %.1f s  %s\n', n, tDesc, files(5).name);
t = tic; readStack(fn, 1:10:n); sec = toc(t);
fprintf('readStack, 1 page in 10: %.1f s for %d of %d pages  %s\n', sec, numel(1:10:n), n, files(5).name);

prev = [results(strcmp({results.method}, 'previous')).MBps];
new = [results(strcmp({results.method}, 'readStack')).MBps];
fprintf('Mean throughput: previous %.0f MB/s, readStack %.0f MB/s (%.1fx)\n', mean(prev), mean(new), mean(new) / mean(prev));
end
