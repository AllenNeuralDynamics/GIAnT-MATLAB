function o = frameOnlineOffsets(off, times, dt)
%FRAMEONLINEOFFSETS The online DMD offset that each frame mostly saw.
%   o = frameOnlineOffsets(off, times, dt) returns, for each line time t in
%   times, the online offsets [x y] that were on the DMD during the frame,
%   averaged with getImage's kernel (weights exp(-|t - l|/dt) for the lines l
%   with |t - l| <= max(3*dt, linesPerCycle)) and rounded: numel(times) x 2.
%   off is getLineOnlineOffsets' output; y is its DMD-timed yEff.
%
%   A frame formed with slap2Frames(..., CompensateOnline=true,
%   RefOffset=o(f,:)) is in the DMD coordinates of the uncompensated frame,
%   but its lines all agree on one offset: a frame without an x/y switch is
%   unchanged, and in a frame that spans one, the lines on the minority offset
%   are moved to the majority's.

x = double(off.x(:));
y = double(off.y(:));
if isfield(off, 'yEff'), y = double(off.yEff(:)); end
n = numel(x);
P = double(off.linesPerCycle);
dtRead = max(3 * dt, P);
o = zeros(numel(times), 2);
for f = 1:numel(times)
    lns = max(1, ceil(times(f) - dtRead)):min(n, floor(times(f) + dtRead));
    w = exp(-abs(times(f) - lns(:)) / dt);
    o(f, :) = round([w' * x(lns), w' * y(lns)] / sum(w));
end
end
