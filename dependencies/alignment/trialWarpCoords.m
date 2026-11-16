function [rT, cT] = trialWarpCoords(warp, t, r, c)
%TRIALWARPCOORDS Trial-image coordinates of common-frame points.
%
%   [RT, CT] = TRIALWARPCOORDS(WARP, T, R, C) maps common-frame points R, C
%   (rows, cols; any equal size) to trial T's registered image. WARP is the
%   across-trial alignment from SILo: with patch centres (WARP.nodes, from
%   alignTrialsPatchwise) the point x maps to x + u_T(x), the patchwise field;
%   without, it maps through the affine map WARP.affine(:,:,T) (2 x 3, rigid
%   or affine), [rT; cT] = A * [r; c; 1].

sz = size(r);
r = r(:); c = c(:);
if isempty(warp.nodes)
    A = warp.affine(:,:,t);
    rT = reshape(A(1,1)*r + A(1,2)*c + A(1,3), sz);
    cT = reshape(A(2,1)*r + A(2,2)*c + A(2,3), sz);
    return
end
d = warp.disp(:,:,t)';
rT = r; cT = c;
chunk = 50000; %points per block, to bound memory for whole images
for i0 = 1:chunk:numel(r)
    ix = i0:min(i0+chunk-1, numel(r));
    u = patchFieldApply(warp.nodes, warp.sigma, r(ix), c(ix), d);
    rT(ix) = r(ix) + u(:,1);
    cT(ix) = c(ix) + u(:,2);
end
rT = reshape(rT, sz);
cT = reshape(cT, sz);
end
