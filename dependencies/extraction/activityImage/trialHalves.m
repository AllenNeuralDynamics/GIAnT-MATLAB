function half = trialHalves(nFrames, blockLen)
%TRIALHALVES Split a trial's frames into two independent halves of equal size.
%
%   HALF = TRIALHALVES(NFRAMES, BLOCKLEN) labels each frame 1 or 2: alternating
%   blocks of BLOCKLEN frames (blocks rather than single frames, because
%   neighbouring downsampled frames share noise). A short or partial last block
%   can leave one half with far fewer frames, so the larger half's last frames
%   are moved to the other half until the two are equal (to a frame). Their
%   mean is then the trial's mean, and the noise of that mean is a quarter of
%   the variance of their difference, as alignTrialsPatchwise assumes.

half = mod(floor((0:nFrames-1)/max(1, blockLen)), 2) + 1;
n1 = nnz(half == 1);
nMove = floor(abs(2*n1 - nFrames)/2);
if nMove > 0
    big = 1 + (2*n1 < nFrames);
    half(find(half == big, nMove, 'last')) = 3 - big;
end
end
