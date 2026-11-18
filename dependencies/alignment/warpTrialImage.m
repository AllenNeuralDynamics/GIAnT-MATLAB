function out = warpTrialImage(IM, warp, t, method)
%WARPTRIALIMAGE Resample a trial image into the common frame.
%
%   OUT = WARPTRIALIMAGE(IM, WARP, T, METHOD) is like warpImageAffine, with
%   OUT(x, :) = IM(trialWarpCoords(WARP, T, x), :) on the same grid as IM
%   (H x W x C): the patchwise field of trial T, or its affine map. Points
%   outside IM or next to NaN are NaN. METHOD is 'linear' (default) or
%   'cubic', falling back to linear where cubic would reach a NaN.

if nargin < 4 || isempty(method)
    method = 'linear';
end
if isempty(warp.nodes)
    out = warpImageAffine(IM, warp.affine(:,:,t), method);
    return
end
[H, W, C] = size(IM, [1 2 3]);
[cc, rr] = meshgrid(1:W, 1:H);
[rT, cT] = trialWarpCoords(warp, t, rr, cc);
out = nan(H, W, C, 'like', IM);
for k = 1:C
    lin = interp2(IM(:,:,k), cT, rT, 'linear');
    if strcmpi(method, 'cubic')
        cub = interp2(IM(:,:,k), cT, rT, 'cubic');
        cub(isnan(cub)) = lin(isnan(cub));
        out(:,:,k) = cub;
    else
        out(:,:,k) = lin;
    end
end
end
