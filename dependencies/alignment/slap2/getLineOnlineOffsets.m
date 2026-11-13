function off = getLineOnlineOffsets(S2)
%GETLINEONLINEOFFSETS Online motion-correction state of every line of a SLAP2 recording.
%   off = getLineOnlineOffsets(S2) reads the line headers of all data files
%   of S2 in one vectorized pass. S2 is a slap2.Slap2DataFile, its
%   hMultiDataFiles, or a single slap2.util.DataFile. The output struct has:
%     x, y           integer DMD offsets (xOffset_pix, yOffset_pix), numLines x 1 int16
%     z              zOffset_um, numLines x 1 single (zeros if the files don't record it)
%     syncPending    logical, bit 5 of the line flags
%     dmdSyncPin     logical, bit 6 of the line flags
%     xySwitchLines  lines whose x or y offset differs from the previous line's
%     medianXY       [x y] median offset over the recording, rounded (the
%                    default reference of slap2Frames' CompensateOnline)
%     linesPerCycle, numLines
%
%   The online correction changes x/y on the first line acquired after the
%   host has sent the new offset to the DMD, which can be any line of a cycle.
%   z steps on every update but is ignored by slap2Frames (the z actuator is
%   slow and the steps are small).
%
%   Flag bits follow the writer (SLAP2FileFormat.h: syncPending is bit 5,
%   dmdSyncPin is bit 6). parseLineHeader.m and getMotionOffsets.m in the
%   slap2 reader label these two bits the other way round.

if isprop(S2, 'hMultiDataFiles') || isfield(S2, 'hMultiDataFiles')
    files = S2.hMultiDataFiles.hDataFiles;
elseif isprop(S2, 'hDataFiles') || isfield(S2, 'hDataFiles')
    files = S2.hDataFiles;
elseif isprop(S2, 'hDataFile') || isfield(S2, 'hDataFile')
    files = S2.hDataFile;
else
    files = S2;
end

parts = cell(numel(files), 5);
P = numel(files(1).lineHeaderIdxs);
for k = 1:numel(files)
    D = files(k);
    assert(numel(D.lineHeaderIdxs) == P, 'getLineOnlineOffsets:linesPerCycle', ...
        'All files must have the same number of lines per cycle.');
    [parts{k, :}] = readFile(D);
end

off = struct();
off.x = vertcat(parts{:, 1});
off.y = vertcat(parts{:, 2});
off.z = vertcat(parts{:, 3});
off.syncPending = vertcat(parts{:, 4});
off.dmdSyncPin = vertcat(parts{:, 5});
off.linesPerCycle = P;
off.numLines = numel(off.x);
off.xySwitchLines = find(diff(off.x) ~= 0 | diff(off.y) ~= 0) + 1;
off.medianXY = round([median(double(off.x)) median(double(off.y))]);
end

function [x, y, z, syncPending, dmdSyncPin] = readFile(D)
% header layout (bytes): 16-19 flags, 20-21 xOffset_pix, 22-23 yOffset_pix,
% 32-35 zOffset_um (optional); rawData is read as int16 words
if ~strcmpi(D.rawData.Format, 'int16')
    D.rawData.Format = 'int16';
end
nC = double(D.numCycles);
cycleOffsets = (0:nC-1) * (double(D.header.bytesPerCycle) / 2);
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
