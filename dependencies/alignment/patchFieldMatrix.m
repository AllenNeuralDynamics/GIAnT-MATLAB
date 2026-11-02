function [L, E, Lplane] = patchFieldMatrix(nodes, sigma, r, c)
%PATCHFIELDMATRIX Displacements at points from displacements at patch centres.
%
%   [L, E, LPLANE] = PATCHFIELDMATRIX(NODES, SIGMA, R, C) returns the
%   numel(R) x K matrix L with u(R, C) = L * d for the displacements d
%   (K x 1 or K x 2) at the K patch centres NODES (K x 2, [row col]), see
%   alignTrialsPatchwise: the plane least-squares fit through the centre
%   values, plus a blend of each centre's departure from that plane with
%   Gaussian weights of SD SIGMA around the centres, normalized to sum to 1.
%   E (K x K) gives the departures, E * d; LPLANE (numel(R) x K) the plane
%   alone. With fewer than 3 centres the plane is their mean.

K = size(nodes, 1);
r = r(:); c = c(:);
if K >= 3
    mu = mean(nodes, 1);
    sc = max(std(nodes, 0, 1)) + eps; %centred, scaled coordinates for conditioning
    P = [(nodes - mu)/sc, ones(K,1)];
    Pq = [(r - mu(1))/sc, (c - mu(2))/sc, ones(numel(r),1)];
else
    P = ones(K, 1);
    Pq = ones(numel(r), 1);
end
Hp = pinv(P); %plane coefficients from centre values
E = eye(K) - P*Hp;
d2 = (r - nodes(:,1)').^2 + (c - nodes(:,2)').^2;
phi = exp(-(d2 - min(d2, [], 2)) / (2*sigma^2)); %relative to the nearest centre, so far points do not underflow
phi = phi ./ sum(phi, 2);
Lplane = Pq*Hp;
L = Lplane + phi*E;
end
