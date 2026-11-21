function [score, nEvents] = sourceKinetics(E, params)
%SOURCEKINETICS Indicator kinetics of each extracted source's transients.
%   [SCORE, NEVENTS] = SOURCEKINETICS(E, PARAMS) averages each source's
%   strongest transients in its least-squares trace (activity channel,
%   E{t}.dF.ls(:,:,1), pooled over the trials in the cell array E) and
%   compares the average one decay time constant after the peak with the
%   indicator's decay:
%
%       SCORE = (avg(L) / avg(0)) / exp(-L / tau),   L = round(tau) frames,
%
%   with tau = PARAMS.tau_s * PARAMS.analyzeHz frames. Transients with the
%   indicator's kinetics give about 1; transients much faster than the
%   indicator (one or two frames) give about 0.
%
%   Each trial's trace is high-passed (1 s moving mean) and filtered with the
%   indicator's decay (an anti-causal exponential, the matched filter for a
%   step rise and exponential decay). Events are maxima of the filtered trace
%   above 5 robust SDs (median absolute deviation over all trials) and at
%   least 2 tau apart, each aligned on the largest high-passed sample within
%   0.6 tau of the maximum. Frames that are NaN (discarded) count as 0 in
%   the filtering, and events whose peak or decay frame was discarded are
%   skipped.
%
%   SCORE is NaN for sources with fewer than 5 events, and for every source
%   when the decay spans fewer than 2 frames (the kinetics are not resolved).
%   NEVENTS is the number of events averaged for each source.

if ~iscell(E)
    E = {E};
end
nSources = 0;
for t = 1:numel(E)
    if hasTraces(E{t})
        nSources = size(E{t}.dF.ls, 1);
        break
    end
end
score = nan(nSources, 1);
nEvents = zeros(nSources, 1);
tauF = params.tau_s * params.analyzeHz;
if nSources == 0 || tauF < 2
    return
end

thresh = 5;                       %event threshold, robust SDs of the filtered trace
minEvents = 5;
L = round(tauF);                  %lag at which the decay is compared
sep = max(1, round(2 * tauF));    %minimum frames between events
w = max(1, round(0.6 * tauF));    %alignment window, frames either side
hpWin = max(3, round(params.analyzeHz)); %high-pass window, frames (1 s)
a = exp(-1 / tauF);

%high-pass and filter each trial
X = {}; Y = {}; V = {};
for t = 1:numel(E)
    if ~hasTraces(E{t})
        continue
    end
    x = double(E{t}.dF.ls(:, :, 1));
    valid = ~isnan(x);
    x = x - movmean(x, hpWin, 2, 'omitnan');
    x(~valid) = 0;
    y = flip(filter(1 - a, [1 -a], flip(x, 2), [], 2), 2);
    X{end+1} = x; %#ok<AGROW>
    Y{end+1} = y; %#ok<AGROW>
    V{end+1} = valid; %#ok<AGROW>
end

%robust SD of each source's filtered trace over all valid frames
sd = nan(nSources, 1);
for s = 1:nSources
    v = cell2mat(cellfun(@(y, valid) y(s, valid(s, :)), Y, V, 'UniformOutput', false));
    if numel(v) > 10
        sd(s) = 1.4826 * median(abs(v - median(v)));
    end
end

sum0 = zeros(nSources, 1);
sumL = zeros(nSources, 1);
for k = 1:numel(Y)
    x = X{k}; y = Y{k}; valid = V{k};
    T = size(y, 2);
    yd = y;
    yd(~valid) = 0;
    isPk = islocalmax(yd, 2, 'MinSeparation', sep) & yd > thresh * sd;
    for s = find(any(isPk, 2))'
        for p = find(isPk(s, :))
            lo = p - w; hi = p + w;
            if lo < 1 || hi > T
                continue
            end
            [~, j] = max(x(s, lo:hi));
            q = lo + j - 1;
            if q + L > T || ~valid(s, q) || ~valid(s, q + L)
                continue
            end
            sum0(s) = sum0(s) + x(s, q);
            sumL(s) = sumL(s) + x(s, q + L);
            nEvents(s) = nEvents(s) + 1;
        end
    end
end
ok = nEvents >= minEvents & sum0 > 0;
score(ok) = (sumL(ok) ./ sum0(ok)) ./ exp(-L / tauF);
end

function tf = hasTraces(Et)
tf = ~isempty(Et) && isstruct(Et) && isfield(Et, 'dF') && isfield(Et.dF, 'ls') ...
    && ~isempty(Et.dF.ls);
end
