function [r, c] = trialWarpInverse(warp, t, rT, cT)
%TRIALWARPINVERSE Common-frame coordinates of points in a trial's image.
%
%   [R, C] = TRIALWARPINVERSE(WARP, T, RT, CT) inverts trialWarpCoords: the
%   common-frame points that trial T's image points RT, CT (rows, cols; any
%   equal size) are read from. Starts from the inverse of the affine part
%   WARP.affine(:,:,T) and, for a patchwise field, iterates x <- x + (y - map(x)),
%   which converges because the field bends by far less than a pixel per pixel.

sz = size(rT);
y = [rT(:) cT(:)]';
A = warp.affine(:,:,t);
x = A(:,1:2) \ (y - A(:,3));
if ~isempty(warp.nodes)
    for it = 1:20
        [rM, cM] = trialWarpCoords(warp, t, x(1,:), x(2,:));
        e = y - [rM(:) cM(:)]';
        x = x + e;
        if max(abs(e), [], 'all') < 1e-6
            break
        end
    end
end
r = reshape(x(1,:), sz);
c = reshape(x(2,:), sz);
end
