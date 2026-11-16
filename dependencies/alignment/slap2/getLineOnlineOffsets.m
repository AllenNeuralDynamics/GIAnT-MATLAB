function off = getLineOnlineOffsets(S2, opts)
%GETLINEONLINEOFFSETS Online motion-correction state of every line of a SLAP2 recording.
%   off = getLineOnlineOffsets(S2) reads the line headers of all data files
%   of S2 in one vectorized pass. S2 is a slap2.Slap2DataFile, its
%   hMultiDataFiles, or a single slap2.util.DataFile. The output struct has:
%     x, y           integer DMD offsets (xOffset_pix, yOffset_pix), numLines x 1 int16
%     z              zOffset_um, numLines x 1 single (zeros if the files don't record it)
%     syncPending    logical, bit 5 of the line flags
%     dmdSyncPin     logical, bit 6 of the line flags
%     yEff           the y offset on the DMD while each line was acquired
%                    (see below), numLines x 1 int16
%     xySwitchLines  lines whose x or y offset differs from the previous line's
%     medianXY       [x y] median offset over the recording, rounded (the
%                    default reference of slap2Frames' CompensateOnline)
%     linesPerCycle, numLines
%
%   The online correction changes x/y in the header on the first line acquired
%   after the host has sent the new offset to the DMD, which can be any line of
%   a cycle. The host sends y first, then x, then toggles the sync pin that
%   starts the new header values. A new x therefore reaches the DMD on the
%   header switch line itself, but a new y reaches it earlier: YLead(1) lines
%   before the switch line when x doesn't change, and YLead(2) lines before
%   it when x changes too (default [2 6], measured on 2025-03 and 2026-10
%   recordings; the line before the first of these is about two-thirds
%   switched and keeps the old y). yEff applies each new y from then on, but
%   not before the previous x/y switch line.
%
%   off = getLineOnlineOffsets(S2, YLead=[a b]) sets these lead times.
%
%   off = getLineOnlineOffsets(S2, Lines=[first last]) reads only the cycles
%   covering lines first-max(YLead)-1 to last+max(YLead)+1, so that a trial of
%   a long continuous recording doesn't read the headers of the whole
%   recording. The output has the same fields and length (all lines of the
%   recording, so that it is indexed by line as before), and is exact on
%   lines first to last; lines outside the cycles read repeat the first or
%   last value read. medianXY is then the median over the lines read, and
%   off.linesRead gives their first and last line.
%   z steps on every update but is ignored by slap2Frames (the z actuator is
%   slow and the steps are small).
%
%   Flag bits follow the writer (SLAP2FileFormat.h: syncPending is bit 5,
%   dmdSyncPin is bit 6). parseLineHeader.m and getMotionOffsets.m in the
%   slap2 reader label these two bits the other way round.

arguments
    S2
    opts.YLead (1,2) double {mustBeNonnegative, mustBeInteger} = [2 6]
    opts.Lines (1,2) double = [1 Inf]
end

if isprop(S2, 'hMultiDataFiles') || isfield(S2, 'hMultiDataFiles')
    files = S2.hMultiDataFiles.hDataFiles;
elseif isprop(S2, 'hDataFiles') || isfield(S2, 'hDataFiles')
    files = S2.hDataFiles;
elseif isprop(S2, 'hDataFile') || isfield(S2, 'hDataFile')
    files = S2.hDataFile;
else
    files = S2;
end

P = numel(files(1).lineHeaderIdxs);
nCycles = arrayfun(@(D) double(D.numCycles), files(:)');
cyclesBefore = [0 cumsum(nCycles(1:end-1))];
numLines = sum(nCycles) * P;
% the cycles to read: those covering the requested lines, widened by the y lead
pad = max(opts.YLead) + 1;
firstCycle = max(1, floor((opts.Lines(1) - pad - 1) / P) + 1);
lastCycle = min(sum(nCycles), ceil((opts.Lines(2) + pad) / P));
assert(firstCycle <= lastCycle, 'getLineOnlineOffsets:lines', ...
    'Lines [%g %g] are outside the recording (%d lines).', opts.Lines(1), opts.Lines(2), numLines);

parts = cell(numel(files), 5);
for k = 1:numel(files)
    D = files(k);
    assert(numel(D.lineHeaderIdxs) == P, 'getLineOnlineOffsets:linesPerCycle', ...
        'All files must have the same number of lines per cycle.');
    cycles = max(1, firstCycle - cyclesBefore(k)):min(nCycles(k), lastCycle - cyclesBefore(k));
    if ~isempty(cycles)
        [parts{k, :}] = readFile(D, cycles);
    end
end

off = struct();
off.x = vertcat(parts{:, 1});
off.y = vertcat(parts{:, 2});
off.z = vertcat(parts{:, 3});
off.syncPending = vertcat(parts{:, 4});
off.dmdSyncPin = vertcat(parts{:, 5});
off.medianXY = round([median(double(off.x)) median(double(off.y))]);
off.linesRead = [(firstCycle - 1) * P + 1, lastCycle * P];
if off.linesRead(1) > 1 || off.linesRead(2) < numLines
    % same length as the recording, so that callers index by line as before
    nBefore = off.linesRead(1) - 1; nAfter = numLines - off.linesRead(2);
    for f = {'x', 'y', 'z', 'syncPending', 'dmdSyncPin'}
        v = off.(f{1});
        off.(f{1}) = [repmat(v(1), nBefore, 1); v; repmat(v(end), nAfter, 1)];
    end
end
off.linesPerCycle = P;
off.numLines = numel(off.x);
off.xySwitchLines = find(diff(off.x) ~= 0 | diff(off.y) ~= 0) + 1;
off.yEff = dmdY(off.x, off.y, off.xySwitchLines, opts.YLead);
off.yLead = opts.YLead;
end

function [x, y, z, syncPending, dmdSyncPin] = readFile(D, cycles)
% header layout (bytes): 16-19 flags, 20-21 xOffset_pix, 22-23 yOffset_pix,
% 32-35 zOffset_um (optional); rawData is read as int16 words
if ~strcmpi(D.rawData.Format, 'int16')
    D.rawData.Format = 'int16';
end
cycleOffsets = (double(cycles(:)') - 1) * (double(D.header.bytesPerCycle) / 2);
idxs = double(D.lineHeaderIdxs(:)) + cycleOffsets;   % word index of each line header, lines x cycles
idxs = idxs(:)';
w = [8; 9; 10; 11] + idxs;
data = reshape(D.rawData.Data(w(:)), 4, []);
flags = typecast(reshape(data(1:2, :), [], 1), 'uint32');
x = data(3, :)';
y = data(4, :)';
if double(D.header.lineHeaderSizeBytes) > 32
    zw = [16; 17] + idxs;
    z = typecast(reshape(D.rawData.Data(zw(:)), [], 1), 'single');
else
    z = zeros(size(x), 'single');
end
syncPending = bitget(flags, 6) == 1;   % bit 5
dmdSyncPin = bitget(flags, 7) == 1;    % bit 6
end

function yEff = dmdY(x, y, sw, lead)
% y on the DMD: each change of the header y moved lead(1) lines earlier, or
% lead(2) when x changes on the same line, but not before the previous switch
yEff = y;
ySw = find(diff(y) ~= 0) + 1;
for s = ySw(:)'
    k = lead(1);
    if x(s) ~= x(s-1), k = lead(2); end
    prev = sw(find(sw < s, 1, 'last'));
    if isempty(prev), prev = 1; end
    yEff(max(s - k, prev):s-1) = y(s);
end
end
