function [E, sources, kinetics] = rejectFastSources(E, sources, params)
%REJECTFASTSOURCES Remove sources whose transients are faster than the indicator.
%   [E, SOURCES, KINETICS] = REJECTFASTSOURCES(E, SOURCES, PARAMS) scores
%   each source with sourceKinetics (about 1 for the indicator's kinetics,
%   about 0 for one- or two-frame transients). When PARAMS.minKineticScore
%   is above 0, sources scoring below it are removed from every trial's
%   extraction results in the cell array E (footprints, dF.ls, dF.denoised,
%   dF.events, F0, SNR) and from SOURCES (R, C, V). Sources that cannot be
%   scored (NaN: fewer than 5 transients, or a frame rate that does not
%   resolve tau_s) are kept. The removed sources were fitted together with
%   the kept ones, so the kept sources' traces are those of the joint fit.
%
%   KINETICS has fields score and nEvents for the kept sources, and
%   rejected.R, rejected.C, rejected.score, rejected.nEvents for the
%   removed ones (seed positions in the consensus frame, 1-based).

[score, nEvents] = sourceKinetics(E, params);
nSources = numel(sources.R);
if numel(score) ~= nSources
    %no trial returned traces; nothing to score
    score = nan(nSources, 1);
    nEvents = zeros(nSources, 1);
end
keep = true(nSources, 1);
if isfield(params, 'minKineticScore') && ~isempty(params.minKineticScore) ...
        && params.minKineticScore > 0
    keep = ~(score < params.minKineticScore);
end

kinetics = struct();
kinetics.score = score(keep);
kinetics.nEvents = nEvents(keep);
kinetics.rejected = struct('R', sources.R(~keep), 'C', sources.C(~keep), ...
    'score', score(~keep), 'nEvents', nEvents(~keep));

if all(keep)
    return
end
sources.R = sources.R(keep);
sources.C = sources.C(keep);
sources.V = sources.V(keep);
for t = 1:numel(E)
    E{t} = subsetTrial(E{t}, keep);
end
fprintf('Kinetics test: removed %d of %d sources (score < %g)\n', ...
    nnz(~keep), nSources, params.minKineticScore);
end

function Et = subsetTrial(Et, keep)
if isempty(Et) || ~isstruct(Et)
    return
end
if isfield(Et, 'footprints') && size(Et.footprints, 2) == numel(keep)
    Et.footprints = Et.footprints(:, keep);
end
if isfield(Et, 'dF')
    for f = {'ls', 'denoised', 'events'}
        if isfield(Et.dF, f{1}) && size(Et.dF.(f{1}), 1) == numel(keep)
            Et.dF.(f{1}) = Et.dF.(f{1})(keep, :, :);
        end
    end
end
if isfield(Et, 'F0') && size(Et.F0, 1) == numel(keep)
    Et.F0 = Et.F0(keep, :, :);
end
if isfield(Et, 'SNR') && size(Et.SNR, 1) == numel(keep)
    Et.SNR = Et.SNR(keep, :);
end
end
