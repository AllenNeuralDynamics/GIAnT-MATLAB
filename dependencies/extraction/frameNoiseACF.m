function acf = frameNoiseACF(S2data, frameLines, dt, zIdx, nLags, mask)
% FRAMENOISEACF  correlation of photon noise between neighbouring analysis frames of SLAP2 data
%
% getImages forms each analysis frame from the raw lines within dtRead of the
% frame's time, weighted by exp(-|t_frame - t_line| / dt). A superpixel is
% sampled once or a few times per cycle, so neighbouring frames share raw
% samples and their noise is correlated. This predicts that correlation from
% the weights, assuming independent noise on raw samples, averaged over the
% raster superpixels of fast-Z plane zIdx.
%
% S2data     : slap2.Slap2DataFile
% frameLines : line index of each analysis frame (as passed to getImages)
% dt         : the dt passed to getImages (lines)
% nLags      : number of lags to return
% mask       : (optional) logical image, the size of getImages' raster
%              output, of the pixels to average over (default all)
% acf        : 1 x nLags correlation at lags 1..nLags (frames)

hdr = S2data.hMultiDataFiles.header;
P = double(hdr.linesPerCycle);
dtRead = max(3 * dt, P);
nRaster = double(hdr.dmdPixelsPerRow) * double(hdr.dmdPixelsPerColumn);

%sampling lines (1..P) of each raster superpixel
lines = find(S2data.hMultiDataFiles.lineFastZIdxs == zIdx);
ids = cell(1, numel(lines));
for lix = 1:numel(lines)
    ids{lix} = double(S2data.hMultiDataFiles.lineSuperPixelIDs{lines(lix)}(:)) + 1;
end
nPer = cellfun(@numel, ids);
spId = vertcat(ids{:});
spLine = repelem(lines(:), nPer(:));
keep = spId <= nRaster;
if nargin > 5 && ~isempty(mask)
    maskT = mask.'; %getImages' raster image is the transpose of superpixel order
    keep(keep) = maskT(spId(keep));
end
spId = spId(keep);
spLine = spLine(keep);
assert(~isempty(spId), 'frameNoiseACF:noPixels', 'no raster superpixels to average over');
[~, ~, g] = unique(spId);
offsets = accumarray(g, spLine, [], @(x) {sort(x)'});

%superpixels sharing sampling lines share the same correlation
key = cellfun(@(x) sprintf('%d,', x), offsets, 'UniformOutput', false);
[~, first, pat] = unique(key);
weight = accumarray(pat, 1);

%a stretch of frames from the middle of the trial
nF = min(numel(frameLines), 200);
i0 = floor((numel(frameLines) - nF) / 2) + 1;
tau = double(frameLines(i0:i0 + nF - 1));
tau = tau(:);
f = (1 + nLags):(nF - 2 * nLags); %frames away from the stretch's edges
cycles = floor((tau(1) - dtRead) / P) - 1 : ceil((tau(end) + dtRead) / P) + 1;

acf = zeros(1, nLags);
for p = 1:numel(first)
    t = offsets{first(p)}(:) + cycles * P;
    t = t(:)';
    d = abs(tau - t);
    W = exp(-d / dt) .* (d < dtRead);
    W = W ./ sum(W, 2);
    C = W * W';
    v = diag(C);
    for L = 1:nLags
        acf(L) = acf(L) + weight(p) * mean(C(sub2ind(size(C), f, f + L))' ./ sqrt(v(f) .* v(f + L)));
    end
end
acf = acf / sum(weight);
end
