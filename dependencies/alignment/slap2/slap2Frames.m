function [IM, freshness] = slap2Frames(S2, channels, times, dt, opts)
%SLAP2FRAMES Raster images from SLAP2 line data, optionally undoing online motion correction.
%   [IM, freshness] = slap2Frames(S2, channels, times, dt) returns, for each
%   line time t in times, the raster image S2.getImage(channel, t, dt, 1, 1)
%   returns. IM is rows x cols x numel(channels) x numel(times) single.
%   freshness, the summed sample weights of the first channel (getImage's
%   third output), is rows x cols x numel(times).
%
%   [...] = slap2Frames(..., Name=Value) options:
%     Kernel            'exp' (default) weights lines by exp(-|t - line|/dt) for
%                       |t - line| <= max(3*dt, linesPerCycle), as getImage does.
%                       'cycle' treats times as cycle numbers. Each frame then
%                       averages all lines of that cycle with equal weight, and
%                       dt is unused.
%     CompensateOnline  false (default). true moves each line's samples by the
%                       integer DMD offset that the online motion correction had
%                       applied to that line, relative to RefOffset, so that all
%                       lines share one DMD frame. Only x/y are undone; z steps
%                       are small and the z actuator is slow.
%     Offsets           output of getLineOnlineOffsets(S2), read when empty.
%                       Pass it in when calling repeatedly on one recording.
%     RefOffset         header offset [x y] that is placed without a shift.
%                       Default: the median offset over the recording.
%     Mapping           2x2 matrix taking a header offset difference [dx; dy]
%                       to the [row; col] displacement applied to the samples.
%                       Default: see defaultMapping below.
%     Rows, Cols        image rows and columns to return (default: all).
%     ZIdx              fast-z plane (default 1).
%
%   Only raster superpixels are used, as with getImage(..., spTypeFlag = 1).
%   The sums are formed by slap2FrameSums, which reads the data once per call,
%   so asking for many frames in one call is much faster than one call per
%   frame; a single frame still takes several times less than getImage.

arguments
    S2
    channels (1,:) double {mustBeInteger, mustBePositive}
    times (1,:) double
    dt (1,1) double {mustBePositive} = 1
    opts.Kernel (1,:) char {mustBeMember(opts.Kernel, {'exp', 'cycle'})} = 'exp'
    opts.CompensateOnline (1,1) logical = false
    opts.Offsets = []
    opts.RefOffset double = []
    opts.Mapping double = []
    opts.Rows (1,:) double = []
    opts.Cols (1,:) double = []
    opts.ZIdx (1,1) double = 1
    opts.Reader (1,:) char = 'auto'
end

hdr = S2.hMultiDataFiles.header;
NR = double(hdr.dmdPixelsPerRow);      % image columns
NC = double(hdr.dmdPixelsPerColumn);   % image rows
rows = opts.Rows; if isempty(rows), rows = 1:NC; end
cols = opts.Cols; if isempty(cols), cols = 1:NR; end
nRows = numel(rows); nCols = numel(cols); nOut = nRows * nCols;
nF = numel(times); nCh = numel(channels);

%only the pixels of the window are formed
window = [];
if nRows < NC || nCols < NR
    [cc, rr] = meshgrid(cols, rows);
    window = (rr(:) - 1) * NR + cc(:);
end
map = opts.Mapping;
if isempty(map), map = defaultMapping(); end
[Num, Den, pix] = slap2FrameSums(S2, channels, times, dt, Kernel = opts.Kernel, ...
    CompensateOnline = opts.CompensateOnline, Offsets = opts.Offsets, RefOffset = opts.RefOffset, ...
    Mapping = map, Pixels = window, ZIdx = opts.ZIdx, Reader = opts.Reader);

%place them in the window
rowPos = zeros(NC, 1); rowPos(rows) = 1:nRows;
colPos = zeros(NR, 1); colPos(cols) = 1:nCols;
pixOut = rowPos(floor((pix - 1) / NR) + 1) + nRows * (colPos(mod(pix - 1, NR) + 1) - 1);
IMw = nan(nOut, nF, nCh, 'single'); Fw = zeros(nOut, nF, 'single');
IMw(pixOut, :, :) = single(Num ./ Den);
if ~isempty(pixOut), Fw(pixOut, :) = single(Den(:, :, 1)); end
IM = permute(reshape(IMw, nRows, nCols, nF, nCh), [1 2 4 3]);
freshness = reshape(Fw, nRows, nCols, nF);
end

function map = defaultMapping()
% Header offset difference [dx; dy] -> [row; col] displacement of a line's
% samples: rows follow yOffset_pix and columns follow xOffset_pix, both with a
% positive sign. Fitted on two 2025-03 trials and two 2026-10 chunks: across
% x/y switches, uncompensated image content moves by about -0.5*[dy; dx] (the
% remainder is tissue motion the step was correcting), and the other signed
% permutations fit worse in every file.
map = [0 1; 1 0];
end
