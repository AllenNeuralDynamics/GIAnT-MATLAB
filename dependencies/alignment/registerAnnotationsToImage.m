function [roiEntry, ok, shiftRC, ncc] = registerAnnotationsToImage(roiEntry, IMtarget, minNCC)
%REGISTERANNOTATIONSTOIMAGE Translate ROI annotations into the space of IMtarget.
%
% [roiEntry, ok, shiftRC, ncc] = registerAnnotationsToImage(roiEntry, IMtarget, minNCC)
%
% roiEntry.IM is the image the ROIs were drawn on (H x W or H x W x C) and
% roiEntry.roiData the ROIs from drawROIs. IMtarget is the image the ROIs
% should be expressed in (e.g. the registered mean image of the first valid
% trial). The integer translation between the two is found by masked normalized
% cross-correlation, and ROI coordinates/masks are shifted so that
%   targetCoord = annotationCoord + shiftRC   ([row col]).
% If the peak correlation is below minNCC (default 0.6), ok is false and the
% ROIs are returned unchanged.

if nargin < 3 || isempty(minNCC)
    minNCC = 0.6;
end
ok = false;
shiftRC = [nan nan];
ncc = nan;

IMann = collapseImage(getFieldOrEmpty(roiEntry, 'IM'));
IMtarget = collapseImage(IMtarget);
if isempty(IMann) || isempty(IMtarget) || ~any(isfinite(IMann(:))) || ~any(isfinite(IMtarget(:)))
    return
end

[shiftRC, ncc] = estimateShift(IMtarget, IMann);
if ~(ncc >= minNCC)
    return
end
ok = true;

targetSz = size(IMtarget, [1 2]);
roiData = getFieldOrEmpty(roiEntry, 'roiData');
for rix = 1:numel(roiData)
    S = roiData{rix};
    for fnm = {'Center', 'Position', 'Vertices'}
        if isfield(S, fnm{1}) && ~isempty(S.(fnm{1}))
            S.(fnm{1}) = S.(fnm{1}) + [shiftRC(2) shiftRC(1)]; %images.roi coordinates are [x y]
        end
    end
    if isfield(S, 'mask') && ~isempty(S.mask)
        S.mask = shiftMask(S.mask, shiftRC, targetSz);
    end
    roiData{rix} = S;
end
roiEntry.roiData = roiData;
end


function [shiftRC, ncc] = estimateShift(IMtarget, IMann)
%offset of IMann's origin within IMtarget by masked normalized cross-correlation
%(Padfield 2012): statistics are computed only over pixels finite in both
%images, so nan regions and partial overlap do not bias the peak. Shifts where
%the overlap is small are excluded, as the correlation is unreliable there.
mT = isfinite(IMtarget); T = IMtarget; T(~mT) = 0;
mM = isfinite(IMann);    M = IMann;    M(~mM) = 0;
outSz = size(T) + size(M) - 1;
xc = @(A, B) real(ifft2(fft2(A, outSz(1), outSz(2)) .* fft2(rot90(B, 2), outSz(1), outSz(2))));

N = round(xc(double(mT), double(mM)));
sumT = xc(T, double(mM));    sumM = xc(double(mT), M);
sumT2 = xc(T.^2, double(mM)); sumM2 = xc(double(mT), M.^2);
sumTM = xc(T, M);
varT = sumT2 - sumT.^2 ./ N;
varM = sumM2 - sumM.^2 ./ N;
C = (sumTM - sumT .* sumM ./ N) ./ sqrt(varT .* varM);

minOverlap = 0.3 * min(nnz(mT), nnz(mM));
C(N < minOverlap | varT <= eps | varM <= eps | ~isfinite(C)) = -inf;

[ncc, ix] = max(C(:));
[r, c] = ind2sub(outSz, ix);
shiftRC = [r - size(M,1), c - size(M,2)];
end


function IM = collapseImage(IMraw)
IM = double(IMraw);
if ndims(IM) > 2
    IM = mean(IM, 3:ndims(IM), 'omitnan');
end
end


function maskOut = shiftMask(maskIn, shiftRC, targetSz)
maskOut = false(targetSz);
[rr, cc] = find(maskIn);
rr = rr + shiftRC(1);
cc = cc + shiftRC(2);
valid = rr >= 1 & rr <= targetSz(1) & cc >= 1 & cc <= targetSz(2);
maskOut(sub2ind(targetSz, rr(valid), cc(valid))) = true;
end


function val = getFieldOrEmpty(S, fieldName)
if isfield(S, fieldName)
    val = S.(fieldName);
else
    val = [];
end
end
