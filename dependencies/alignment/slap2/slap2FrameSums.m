function [Num, Den, pix] = slap2FrameSums(S2, channels, times, dt, opts)
%SLAP2FRAMESUMS Weighted sums of SLAP2 raster samples for a set of frames.
%   [Num, Den, pix] = slap2FrameSums(S2, channels, times, dt) forms, for each
%   line time t in times, the sums that getImage and getImages divide, at the
%   raster superpixels that the frames' lines sample:
%     Num(i, f, k)  sum over lines l of w(f, l) * (sample of channel channels(k)
%                   at raster superpixel pix(i) on line l)
%     Den(i, f, k)  sum of w(f, l) over the same samples
%   so that Num ./ Den is the image value. pix (nPix x 1, sorted) are 1-based
%   raster superpixel IDs: image column mod(pix-1, NR)+1 and row
%   floor((pix-1)/NR)+1, with NR = dmdPixelsPerRow. slap2Frames forms images
%   from these sums, and rawFramesAt (extraction) values at chosen pixels.
%
%   [...] = slap2FrameSums(..., Name=Value) options:
%     Kernel            'exp' (default): w(f, l) = exp(-|t - l|/dt) for lines with
%                       |t - l| <= dtRead = max(3*dt, linesPerCycle), as getImage.
%                       'cycle': times are cycle numbers, and each frame weights
%                       all lines of its cycle equally (dt is unused).
%     Boundary          'inclusive' (default) as getImage; 'exclusive' leaves out
%                       the lines exactly dtRead from t, as getImages does.
%     MaskInvalid       true (default) skips samples <= intmin('int16'), as
%                       getImage does; Den is then per channel. false sums them,
%                       as getImages does, and Den is nPix x nF.
%     CompensateOnline  false (default). true moves each line's samples by the
%                       integer DMD offset that the online motion correction had
%                       applied to that line, relative to RefOffset, so that all
%                       lines share one DMD frame; pix are then positions in that
%                       frame. Only x/y are undone: z steps are small and the z
%                       actuator is slow.
%     Offsets           output of getLineOnlineOffsets(S2), read when empty.
%     RefOffset         header offset [x y] placed without a shift. Default:
%                       the median offset over the recording, rounded.
%     Mapping           2x2 matrix taking a header offset difference [dx; dy] to
%                       the [row; col] displacement of the samples. Default
%                       [0 1; 1 0]: rows follow y and columns follow x (measured;
%                       see slap2Frames).
%     Pixels            raster IDs to form (default: all that are sampled); pix
%                       is then the subset of these that is sampled.
%     ZIdx              fast-z plane (default 1).
%
%   Each line index samples the same superpixels in every cycle, so for each
%   line index the samples of all cycles in the request form a dense matrix
%   (superpixels x lines). One product with that line index's kernel weights
%   (lines x frames) gives its contribution to every frame, and a sparse slot ->
%   pixel map (one per online offset present) adds the results up per pixel.
%   The cycles are read from disk in one contiguous block per data file, which
%   is several times faster than the reader's getLineData; objects without data
%   files (test doubles) are read through getLineData.

arguments
    S2
    channels (1,:) double {mustBeInteger, mustBePositive}
    times (1,:) double
    dt (1,1) double {mustBePositive} = 1
    opts.Kernel (1,:) char {mustBeMember(opts.Kernel, {'exp', 'cycle'})} = 'exp'
    opts.Boundary (1,:) char {mustBeMember(opts.Boundary, {'inclusive', 'exclusive'})} = 'inclusive'
    opts.MaskInvalid (1,1) logical = true
    opts.CompensateOnline (1,1) logical = false
    opts.Offsets = []
    opts.RefOffset double = []
    opts.Mapping double = []
    opts.Pixels double = []
    opts.ZIdx (1,1) double = 1
    opts.Reader (1,:) char {mustBeMember(opts.Reader, {'auto', 'file', 'getLineData'})} = 'auto'
end

M = S2.hMultiDataFiles;
hdr = M.header;
P = double(hdr.linesPerCycle);
NR = double(hdr.dmdPixelsPerRow);      % image columns
NC = double(hdr.dmdPixelsPerColumn);   % image rows
nPx = NR * NC;
totalLines = double(M.totalNumLines);
nF = numel(times); nCh = numel(channels);
nDen = 1; if opts.MaskInvalid, nDen = nCh; end

%lines that contribute to each frame, and their weights (as getImage / getImages)
frameLines = cell(1, nF); frameWeights = cell(1, nF);
for f = 1:nF
    switch opts.Kernel
        case 'exp'
            dtRead = max(3 * dt, P);
            lns = max(1, times(f) - dtRead):min(totalLines, times(f) + dtRead);
            if strcmp(opts.Boundary, 'exclusive')
                lns = lns(abs(times(f) - lns) < dtRead);
            end
            wts = exp(-abs((times(f) - lns) / dt));
        case 'cycle'
            lns = ((times(f) - 1) * P + 1):min(totalLines, times(f) * P);
            wts = ones(size(lns));
    end
    assert(all(lns == round(lns)), 'slap2Frames:nonIntegerLines', ...
        'times and dt must give integer line numbers.');
    keep = M.lineFastZIdxs(mod(lns - 1, P) + 1) == opts.ZIdx;
    frameLines{f} = lns(keep(:)'); frameWeights{f} = wts(keep(:)');
end
uLines = unique([frameLines{:}]);
nU = numel(uLines);
if nU == 0
    Num = zeros(0, nF, nCh); Den = zeros(0, nF, nDen); pix = zeros(0, 1);
    return
end
lookup = zeros(1, uLines(end) - uLines(1) + 1);
lookup(uLines - uLines(1) + 1) = 1:nU;
wi = cell(1, nF); wj = cell(1, nF);
for f = 1:nF
    wi{f} = lookup(frameLines{f} - uLines(1) + 1);
    wj{f} = f * ones(1, numel(frameLines{f}));
end
Wt = sparse([wj{:}], [wi{:}], [frameWeights{:}], nF, nU);   % frames x lines
li = mod(uLines(:) - 1, P) + 1; ci = floor((uLines(:) - 1) / P) + 1;

tbl = slotTable(M, P, NR, nPx);
nSlots = tbl.nSlots; nSlotL = tbl.nSlotL; slotStart = tbl.slotStart;

%integer shift of each line's samples, as an index into the distinct shifts
shifts = [0 0]; shiftIdx = ones(nU, 1);
if opts.CompensateOnline
    off = opts.Offsets;
    if isempty(off), off = getLineOnlineOffsets(S2); end
    ref = opts.RefOffset;
    if isempty(ref) && isfield(off, 'medianXY'), ref = off.medianXY; end
    if isempty(ref), ref = round([median(double(off.x)) median(double(off.y))]); end
    map = opts.Mapping;
    if isempty(map), map = [0 1; 1 0]; end
    d = (map * [double(off.x(uLines))' - ref(1); double(off.y(uLines))' - ref(2)])';
    assert(all(d(:) == round(d(:))), 'slap2Frames:nonIntegerShift', 'Mapping must give integer shifts.');
    [shifts, ~, shiftIdx] = unique(d, 'rows');
end
nS = size(shifts, 1);

%where each shift's sums go. Without a pixel subset, the sums of the lines of
%each shift are formed at the unshifted pixels pix0 with the cached slot maps,
%and then moved by the shift: dest(i, s) is the output row of pix0(i) under
%shift s (0: off the image). With a subset, only the slots that land on wanted
%pixels are formed: Q(slot, s) is the output row of each slot under shift s.
translate = isempty(opts.Pixels);
if translate
    pr = floor((tbl.pix0 - 1) / NR) + 1; pc = mod(tbl.pix0 - 1, NR) + 1;
    sid = zeros(numel(tbl.pix0), nS);
    for s = 1:nS
        r = pr + shifts(s, 1); c = pc + shifts(s, 2);
        ok = r >= 1 & r <= NC & c >= 1 & c <= NR;
        sid(ok, s) = (r(ok) - 1) * NR + c(ok);
    end
    hit = sid > 0;
    dest = zeros(size(sid));
    if nS == 1   % a translation keeps the pixels sorted
        pix = sid(hit); dest(hit) = 1:numel(pix);
    else
        mark = false(nPx, 1); mark(sid(hit)) = true;
        pix = find(mark);
        lut = zeros(nPx, 1); lut(pix) = 1:numel(pix);
        dest(hit) = lut(sid(hit));
    end
else
    want = false(nPx, 1); want(opts.Pixels) = true;
    q = zeros(nSlots, nS);
    for s = 1:nS
        r = tbl.r0 + shifts(s, 1); c = tbl.c0 + shifts(s, 2);
        ok = r >= 1 & r <= NC & c >= 1 & c <= NR;
        id = (r(ok) - 1) * NR + c(ok);
        good = want(id);
        okIdx = find(ok);
        q(okIdx(good), s) = id(good);
    end
    hit = q > 0;
    mark = false(nPx, 1); mark(q(hit)) = true;
    pix = find(mark);
    lut = zeros(nPx, 1); lut(pix) = 1:numel(pix);
    Q = zeros(size(q)); Q(hit) = lut(q(hit));
end
nPix = numel(pix);

%the data: one contiguous block of whole cycles, or the reader's cells
useFile = strcmp(opts.Reader, 'file') || (strcmp(opts.Reader, 'auto') && hasDataFiles(M));
if useFile
    c1 = min(ci);
    block = readCycles(M, c1, max(ci));
    rowBase = double(M.lineDataStartIdxs(:)) - double(M.lineHeaderIdxs(1));
    nPerCh = double(M.lineDataNumElements(:)) / double(hdr.numChannels);
else
    %the reader returns channels min:max in file order
    chRange = min(channels):max(channels);
    chCols = channels - chRange(1) + 1;
    ld = M.getLineData(li, ci, chRange);
end

%per shift and line index: samples (slots x lines) times weights (lines x frames).
%Each shift only touches the frames whose kernels overlap its lines, and the
%denominator is (slots per pixel and line index) x (summed weights), corrected
%only where samples are invalid.
Num = zeros(nPix, nF, nCh); Den = zeros(nPix, nF, nDen);
[~, order] = sortrows([shiftIdx li]);
grpKey = (shiftIdx(order) - 1) * P + li(order);
edges = [find([true; diff(grpKey) ~= 0]); nU + 1];
grpShift = shiftIdx(order(edges(1:end-1)));
invalid = double(intmin('int16'));
for s = 1:nS
    if translate
        rowsOut = find(dest(:, s) > 0);           % rows of the sums below that stay on the image
        if isempty(rowsOut), continue, end
        dst = dest(rowsOut, s);
        nm = nSlots; pos = []; Gs = tbl.G0; Cs = tbl.Cs0;
    else
        mapped = find(Q(:, s) > 0);
        nm = numel(mapped);
        if nm == 0, continue, end
        pos = zeros(nSlots, 1); pos(mapped) = 1:nm;
        Gs = sparse(Q(mapped, s), (1:nm)', 1, nPix, nm);               % slot -> pixel
        Cs = sparse(Q(mapped, s), tbl.slotLine(mapped), 1, nPix, P);   % pixel x line index: slots landing there
        rowsOut = (1:nPix)'; dst = rowsOut;
    end
    fcols = find(any(Wt(:, shiftIdx == s), 2))';
    nFs = numel(fcols);
    Y = zeros(nm, nFs, nCh); Ssum = zeros(P, nFs); Vbad = [];
    for g = find(grpShift == s)'
        members = order(edges(g):edges(g+1) - 1);
        L = li(members(1));
        Wg = full(Wt(fcols, members)).';                               % lines x frames
        Ssum(L, :) = sum(Wg, 1);
        if translate
            rowsY = slotStart(L) + (1:nSlotL(L))';
            sel = tbl.slotSel{L};
        else
            within = find(pos(slotStart(L) + (1:nSlotL(L))) > 0);
            rowsY = pos(slotStart(L) + within);
            sel = tbl.slotSel{L}(within);
        end
        if isempty(sel), continue, end
        if ~useFile
            X = cat(3, ld{members});                                   % superpixels x channels x lines
        end
        for k = 1:nCh
            if useFile
                Dk = double(block(rowBase(L) + (channels(k) - 1) * nPerCh(L) + sel, ci(members) - c1 + 1));
            else
                Dk = double(reshape(X(sel, chCols(k), :), numel(sel), []));
            end
            if opts.MaskInvalid
                bad = Dk <= invalid;
                if any(bad(:))
                    if isempty(Vbad), Vbad = zeros(nm, nFs, nCh); end
                    Dk(bad) = 0;
                    Vbad(rowsY, :, k) = double(bad) * Wg;
                end
            end
            Y(rowsY, :, k) = Dk * Wg;
        end
    end
    for k = 1:nCh
        add = Gs * Y(:, :, k);
        Num(dst, fcols, k) = Num(dst, fcols, k) + add(rowsOut, :);
    end
    dk = Cs * Ssum;
    for k = 1:nDen
        add = dk;
        if ~isempty(Vbad), add = add - Gs * Vbad(:, :, k); end
        Den(dst, fcols, k) = Den(dst, fcols, k) + add(rowsOut, :);
    end
end
end

function tbl = slotTable(M, P, NR, nPx)
% raster slots of the cycle (line index, position in the line's data, image
% row/col) and their maps onto the sampled pixels without shifts; the last
% recording's table is kept, since images are usually formed one frame or block
% at a time from one recording
persistent cache
key = '';
if isprop(M, 'filename') || isfield(M, 'filename')
    key = sprintf('%s|%d|%d', M.filename, P, numel(M.lineSuperPixelIDs));
end
if ~isempty(key) && ~isempty(cache) && strcmp(cache.key, key)
    tbl = cache.tbl;
    return
end
tbl.slotSel = cell(P, 1); slotId = cell(P, 1);
for L = 1:P
    ids = double(M.lineSuperPixelIDs{L}(:)) + 1;
    sel = find(ids <= nPx);
    tbl.slotSel{L} = sel;
    slotId{L} = ids(sel);
end
tbl.nSlotL = cellfun(@numel, tbl.slotSel);
tbl.slotStart = [0; cumsum(tbl.nSlotL)];
tbl.nSlots = tbl.slotStart(end);
baseId = vertcat(zeros(0, 1), slotId{:});
tbl.r0 = floor((baseId - 1) / NR) + 1;
tbl.c0 = mod(baseId - 1, NR) + 1;
tbl.slotLine = repelem((1:P)', tbl.nSlotL);
[tbl.pix0, ~, Q0] = unique(baseId);
tbl.Q0 = Q0(:);
nPix0 = numel(tbl.pix0);
tbl.G0 = sparse(tbl.Q0, (1:tbl.nSlots)', 1, nPix0, tbl.nSlots);
tbl.Cs0 = sparse(tbl.Q0, tbl.slotLine, 1, nPix0, P);
if ~isempty(key)
    cache = struct('key', key, 'tbl', tbl);
end
end

function tf = hasDataFiles(M)
tf = (isprop(M, 'hDataFiles') || isfield(M, 'hDataFiles')) && ~isempty(M.hDataFiles) ...
    && (isprop(M.hDataFiles(1), 'datFileName') || isfield(M.hDataFiles(1), 'datFileName'));
end

function block = readCycles(M, firstCycle, lastCycle)
% the raw int16 words of cycles firstCycle..lastCycle (across data files), one
% column per cycle; word lineDataStartIdxs(L) - lineHeaderIdxs(1) + j of a
% column is word j of line L's data in that cycle
assert(~M.header.channelsInterleave, 'Only deinterleaved data is supported');
files = M.hDataFiles;
fileEnd = cumsum(double([files.numCycles]));
fileStart = [0, fileEnd(1:end-1)] + 1;
parts = {};
for f = find(fileStart <= lastCycle & fileEnd >= firstCycle)
    D = files(f);
    c1 = max(firstCycle, fileStart(f)) - fileStart(f) + 1;
    c2 = min(lastCycle, fileEnd(f)) - fileStart(f) + 1;
    bytesPerCycle = double(D.header.bytesPerCycle);
    fid = fopen(D.datFileName, 'r');
    assert(fid > 0, 'slap2Frames:readFailed', 'Failed to open %s.', D.datFileName);
    closeFile = onCleanup(@() fclose(fid));
    fseek(fid, 2 * (double(D.lineHeaderIdxs(1)) - 1) + (c1 - 1) * bytesPerCycle, 'bof');
    parts{end+1} = fread(fid, [bytesPerCycle / 2, c2 - c1 + 1], '*int16'); %#ok<AGROW>
    assert(size(parts{end}, 2) == c2 - c1 + 1, 'slap2Frames:readFailed', 'Failed to read cycles of %s.', D.datFileName);
    clear closeFile
end
block = [parts{:}];
end
