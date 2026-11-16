function u = patchFieldApply(nodes, sigma, r, c, d)
%PATCHFIELDAPPLY Patchwise field at points: patchFieldMatrix(NODES, SIGMA, R, C) * D.
%
%   U = PATCHFIELDAPPLY(NODES, SIGMA, R, C, D) returns the displacements U
%   (numel(R) x size(D, 2)) at the points R, C of the field set by the values D
%   (K x 1 or K x 2) at the K patch centres NODES, without forming the
%   numel(R) x K field matrix: the plane coefficients and the centres'
%   departures from the plane are computed from D once, then the plane is
%   evaluated and the departures blended at each point. Same result as
%   patchFieldMatrix(NODES, SIGMA, R, C) * D up to rounding.

K = size(nodes, 1);
r = r(:); c = c(:);
if K >= 3
    mu = mean(nodes, 1);
    sc = max(std(nodes, 0, 1)) + eps;
    P = [(nodes - mu)/sc, ones(K,1)];
    Pq = [(r - mu(1))/sc, (c - mu(2))/sc, ones(numel(r),1)];
else
    P = ones(K, 1);
    Pq = ones(numel(r), 1);
end
coef = pinv(P) * d;      % plane coefficients
dep = d - P * coef;      % departures of the centre values from the plane
d2 = (r - nodes(:,1)').^2 + (c - nodes(:,2)').^2;
phi = exp(-(d2 - min(d2, [], 2)) / (2*sigma^2));
phi = phi ./ sum(phi, 2);
u = Pq * coef + phi * dep;
end
