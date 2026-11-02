function [IM, vIM] = interpFramesAt(M1, R, C, freshness)
%INTERPFRAMESAT interpolates a multicolor image at arbitrary points and computes
%the expected variance multiplier vIM
%
%   Like interpFrames, but the query points R, C (rows, cols in M1; vectors of
%   equal length) need not form a separable grid, so they can follow an affine
%   (or any) map. Bilinear, weighted by freshness; a point is NaN if any of its
%   four neighbours is outside M1 or NaN. IM is numel(R) x channels.
sz = size(M1, [1 2]);
nCh = size(M1, 3);
R = R(:); C = C(:);
fr = mod(R,1); fc = mod(C,1);
c1 = (1-fc).*(1-fr);
c2 = fc.*(1-fr);
c3 = (1-fc).*fr;
c4 = fc.*fr;

[a1, f1] = gatherCorner(floor(R), floor(C));
[a2, f2] = gatherCorner(floor(R), ceil(C));
[a3, f3] = gatherCorner(ceil(R), floor(C));
[a4, f4] = gatherCorner(ceil(R), ceil(C));

%the variance of each a1...a4 is (n1*a1)/n1.^2, i.e. (original photon
%counts)/averaging factor
den = c1.*f1 + c2.*f2 + c3.*f3 + c4.*f4;
IM = (c1.*f1.*a1 + c2.*f2.*a2 + c3.*f3.*a3 + c4.*f4.*a4) ./ den;
vIM = 1./den; %this is not the single-pixel variance, but rather a local variance multiplier
vIM(isinf(vIM)) = nan;

    function [a, f] = gatherCorner(rr, cc)
        invalid = rr<1 | rr>sz(1) | cc<1 | cc>sz(2);
        rr(invalid) = 1; cc(invalid) = 1;
        ix = rr + (cc-1)*sz(1);
        a = reshape(M1(ix + (0:nCh-1)*prod(sz)), [], nCh);
        a(invalid,:) = nan;
        f = freshness(ix);
        f(invalid) = nan;
    end
end
