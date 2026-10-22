function roiData = shiftRoiData(roiData, shiftRC, targetSz)
%SHIFTROIDATA Translate drawROIs ROIs by shiftRC = [row col] pixels.
%
% roiData = shiftRoiData(roiData, shiftRC, targetSz)
%
% A point at [row col] moves to [row col] + shiftRC. Geometry (Center,
% Position, Vertices; images.roi [x y] order) is shifted exactly; masks are
% shifted by the rounded offset into a targetSz canvas, and pixels that fall
% outside it are dropped.

maskShift = round(shiftRC);
for rix = 1:numel(roiData)
    S = roiData{rix};
    for fnm = {'Center', 'Position', 'Vertices'}
        if isfield(S, fnm{1}) && ~isempty(S.(fnm{1}))
            S.(fnm{1}) = S.(fnm{1}) + [shiftRC(2) shiftRC(1)];
        end
    end
    if isfield(S, 'mask') && ~isempty(S.mask)
        S.mask = shiftMask(S.mask, maskShift, targetSz);
    end
    roiData{rix} = S;
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
