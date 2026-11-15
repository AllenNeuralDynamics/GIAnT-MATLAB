function out = warpImageAffine(IM, A, method)
%WARPIMAGEAFFINE Resample an image into the common frame through an affine map.
%
%   OUT = WARPIMAGEAFFINE(IM, A, METHOD) returns OUT(row,col,:) = IM(A*[row; col; 1], :)
%   on the same grid as IM (H x W x C). A is 2 x 3, mapping common-frame
%   [row; col; 1] to IM's [row; col]. Points outside IM or next to NaN are NaN.
%   METHOD is 'linear' (default) or 'cubic'; 'cubic' falls back to linear
%   where its larger neighbourhood would reach a NaN, so it removes no more
%   pixels at the edges of the scanned region than 'linear'.

if nargin < 3 || isempty(method)
    method = 'linear';
end
[H, W, C] = size(IM, [1 2 3]);
[cc, rr] = meshgrid(1:W, 1:H);
rT = A(1,1)*rr + A(1,2)*cc + A(1,3);
cT = A(2,1)*rr + A(2,2)*cc + A(2,3);
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
