function nfail = testReadStack(dataDir)
%TESTREADSTACK Check that readStack returns exactly what the direct reader returns.
%
% nfail = testReadStack(dataDir)
%
% dataDir: a folder with ScanImage tiffs (up to two of them, <= 3 GB, are
% used), on whatever storage should be tested (local, SMB, S3/FUSE, NFS).
% Also tests h5 movies, sparse parallel reads, compressed and multi-sample
% tiffs, using files written to tempdir. Prints PASS/FAIL lines; returns
% the failure count.

nfail = 0;
fprintf('storageClass(%s) = %s\n', dataDir, storageClass(fullfile(dataDir, 'x.tif')));

%% real tiffs: every strategy, whole stack and a subset of planes
files = dir(fullfile(dataDir, '**', '*.tif'));
files = files([files.bytes] <= 3e9);
[~, o] = sort([files.bytes], 'descend');
files = files(o(1:min(2, end)));
for k = 1:numel(files)
    fn = fullfile(files(k).folder, files(k).name);
    r = ScanImageTiffReader(fn); ref = r.data(); refMeta = r.metadata(); refDesc = r.descriptions(); r.close();
    n = size(ref, ndims(ref));
    sub = unique([1:7:n n]);
    sub = sub([end 1:end-1]); %out of file order
    for strat = {'direct', 'copy', 'parallel', 'auto'}
        [data, meta, desc, info] = readStack(fn, [], 'Strategy', strat{1});
        ok = isequal(data, ref) && isequal(meta, refMeta) && isequal(desc, refDesc) ...
            && (isempty(info.copy) || ~exist(info.copy, 'file'));
        nfail = nfail + check(ok, sprintf('%s whole, %s (used %s on %s)', files(k).name, strat{1}, info.strategy, info.storage));
        [data, ~, desc, info] = readStack(fn, sub, 'Strategy', strat{1});
        ok = isequal(data, ref(:,:,sub)) && isequal(desc, refDesc(sub)) && (isempty(info.copy) || ~exist(info.copy, 'file'));
        nfail = nfail + check(ok, sprintf('%s %d planes, %s (used %s)', files(k).name, numel(sub), strat{1}, info.strategy));
    end
    ok = isequal(ScanImageTiffWrapper(fn), ref);
    nfail = nfail + check(ok, sprintf('%s via ScanImageTiffWrapper', files(k).name));
    clear ref data
end
if isempty(files)
    fprintf('No tiffs <= 3 GB in %s; skipped real-file tests.\n', dataDir);
end

%% h5 movie: whole, subset, named dataset
h5 = [tempname '.h5'];
cleanupH5 = onCleanup(@() delete(h5));
mov = int16(reshape(1:50*40*30, 50, 40, 30));
h5create(h5, '/mov', size(mov), 'Datatype', 'int16'); h5write(h5, '/mov', mov);
nfail = nfail + check(isequal(readStack(h5), mov), 'h5 whole');
nfail = nfail + check(isequal(readStack(h5, [3 8 9 21]), mov(:,:,[3 8 9 21])), 'h5 subset');
nfail = nfail + check(isequal(readStack(h5, 5:6, 'Dataset', 'mov'), mov(:,:,5:6)), 'h5 named dataset');

%% sparse planes in parallel: more reads than workers
big = uint16(randi(4000, 512, 512, 200));
bigFn = [tempname '.tif'];
cleanupBig = onCleanup(@() delete(bigFn));
for p = 1:200, imwrite(big(:,:,p), bigFn, 'Compression', 'none', 'WriteMode', pick(p == 1, 'overwrite', 'append')); end
sel = 200:-20:1; %gaps over 8 MB, so one read per plane
lastwarn('');
data = readStack(bigFn, sel, 'Strategy', 'parallel', 'NumParallel', 2);
[~, wid] = lastwarn; %a sequential or imread fallback would warn
nfail = nfail + check(isequal(data, permute(big(:,:,sel), [2 1 3])) && isempty(wid), ...
    sprintf('%d sparse planes, 2 parallel reads%s', numel(sel), pick(isempty(wid), '', [' (' wid ')'])));
clear big

%% compressed tiff falls back to imread, same layout, keeping the reader's metadata
pages = uint16(randi(4000, 30, 20, 4));
lzw = [tempname '.tif'];
cleanupLzw = onCleanup(@() delete(lzw));
for p = 1:4, imwrite(pages(:,:,p), lzw, 'Compression', 'lzw', 'WriteMode', pick(p == 1, 'overwrite', 'append')); end
legacyMeta = 'SI.hChannels.channelSave = [1;2]';
fid = fopen(lzw, 'a'); fwrite(fid, legacyMeta, 'char'); fwrite(fid, numel(legacyMeta), 'uint32'); fclose(fid); %legacy appended block
s = warning('off', 'readStack:Fallback'); [data, meta] = readStack(lzw); warning(s);
nfail = nfail + check(isequal(data, permute(pages, [2 1 3])) && isequal(meta, legacyMeta), ...
    'compressed tiff via imread fallback, with metadata');

%% multi-sample (RGB) tiff: reader and fallback agree on layout
rgb = uint8(randi(255, 12, 9, 3, 3));
rgbFn = [tempname '.tif'];
cleanupRgb = onCleanup(@() delete(rgbFn));
for p = 1:3, imwrite(rgb(:,:,:,p), rgbFn, 'Compression', 'none', 'WriteMode', pick(p == 1, 'overwrite', 'append')); end
direct = readStack(rgbFn);
fallback = readTiffStackImread(rgbFn);
nfail = nfail + check(isequal(direct, fallback) && isequal(size(direct), [3 9 12 3]), ...
    sprintf('multi-sample tiff: reader %s == imread fallback %s', mat2str(size(direct)), mat2str(size(fallback))));

fprintf('%d failure(s)\n', nfail);
end


function failed = check(ok, msg)
if ok, fprintf('PASS  %s\n', msg); else, fprintf('FAIL  %s\n', msg); end
failed = ~ok;
end


function v = pick(c, a, b)
if c, v = a; else, v = b; end
end
