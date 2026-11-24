function R = subproblemRefitTest(files, varargin)
%SUBPROBLEMREFITTEST Ground-truth test of source extraction on saved subproblems.
%   R = SUBPROBLEMREFITTEST(files) tests how well source extraction recovers
%   activity in the sampling regime of real data. files are subproblem
%   inputs saved by extractTrial with params.saveSubproblemsDir set: a
%   folder (searched recursively for sub_*.mat), a text file listing them
%   one per line, or a cell or string array of file names. For each one:
%     1. Fit it with extractSubproblem (the production solver, unchanged)
%        and take the fit as the truth: baseline B, footprints H and
%        activity X = convn(S, params.k, 'same') (channel 1, the activity
%        channel; a second channel is left out).
%     2. Draw new data from the truth under the fit's own noise model, with
%        the subproblem's freshness F (F_inv_p):
%            Z = F .* poissrnd(T ./ F),   T = max(B + H*X, 0)
%        so Z has mean T and variance F.*T.
%     3. Fit Z again (once per variant, all on the same Z) and compare the
%        refit with the truth.
%   Amplitudes are in the fit's normalized units (the data divided by
%   params.photonScale), at each source's peak pixel: an event of size a
%   adds a to the peak pixel at the event's peak.
%
%   Events are the peaks of each source's true X (local maxima above 2% of
%   its largest, at least L = ceil(2.5*params.tau_full) samples apart; not
%   within 5 samples of the start or L+3 of the end), each with the window
%   from 2 samples before its peak to L after. An event is found when the
%   refit's activity mass in its window is above 25% of the truth's.
%
%   Name-value options:
%     'Variants'     params fields that override the saved params in the
%                    refit: a struct, a struct array or a cell array of
%                    structs, one per variant (default: one, no overrides).
%                    Fields derived from an overridden one (e.g. k and
%                    tau_full from tau_s) are recomputed; photonScale, dXY
%                    and selRadius are fixed by the saved subproblem
%     'VariantNames' cell array of names for the variants (default refit1, ...)
%     'TruthParams'  params fields that override the saved params in the
%                    truth fit (default none)
%     'TruthScale'   multiplies the true activity, e.g. 0.5 or 0.25 to map
%                    the detection threshold (default 1)
%     'InjectSizes'  sizes of events added to the truth at quiet times (no
%                    true event within 2L+4 samples, and less than 10% of
%                    the event's mass already in its window); default none
%     'InjectRate'   injected events per source and second, the sizes in
%                    turn (default 0.5)
%     'Seed'         subproblem i (in the order of files) is simulated with
%                    rng(Seed + i); the generator's state is restored after
%                    (default 0)
%     'SizeEdges'    edges of the event-size bins (default [0 0.25 0.5 1 2 4 8 Inf])
%     'UseParallel'  fit the subproblems in a parfor loop (default false).
%                    Results differ slightly between the two settings, with
%                    the same Seed: a pool worker computes in a single thread
%                    and the client in several, so the fits (and with them
%                    the simulated data) differ in floating-point order.
%                    Compare variants within one call, which refits them all
%                    on the same data; production solves subproblems on
%                    workers
%     'Verbose'      print the report (default true)
%
%   R has fields
%     files, options, variants (the variant names)
%     sources   table, one row per source and variant: file, source,
%               variant, trueMass, recovered (refit / true activity mass),
%               slope (of the refit on the truth, through 0), F0ratio (mean
%               refit F0 / mean true F0, F0 as extractTrial reports it),
%               F0ratioActive (the same where the truth is above 20% of its
%               largest), lsRecovered (least-squares activity mass / true
%               activity mass), peakH (the footprint's largest value) and
%               baseline (the true baseline under the footprint, per pixel
%               and sample)
%     events    table, one row per event and variant: file, source,
%               variant, size, trueMass, refitMass, found, injected
%     summary   table, one row per variant: nSources, nEvents, the
%               median and quartiles of recovered and slope, pooled
%               (refit / true mass summed over events), and medians of
%               F0ratio, F0ratioActive and lsRecovered
%     bySize    table, one row per variant and size bin with events: n,
%               shareFound, recoveredAll (refit / true mass over the bin's
%               events), recoveredFound (over its found events only);
%               injected events are in rows of their own (injected true)
%
%   Example: compare a refit option with the default on the same data
%     R = subproblemRefitTest('D:\subproblems\path1_trial002_x', ...
%         'Variants', {struct(), struct('nmfIter', 3)}, 'VariantNames', {'default', 'nmfIter3'});
%
%   Limits: the truth uses the fit's own footprints and kernel, so the test
%   cannot show a footprint-width or kinetics mismatch; noise is drawn
%   independently per sample.
%
%   See also extractSubproblem, extractTrial.

opts = parseOptions(varargin{:});
fileList = listFiles(files);
if isempty(fileList)
    error('subproblemRefitTest:noFiles', 'No subproblem files (sub_*.mat) found.');
end
nV = numel(opts.Variants);
res = cell(numel(fileList), 1);
t0 = tic;
if opts.UseParallel
    parfor i = 1:numel(fileList)
        res{i} = oneSubproblem(fileList{i}, i, opts);
    end
else
    for i = 1:numel(fileList)
        res{i} = oneSubproblem(fileList{i}, i, opts);
        if opts.Verbose
            fprintf('%d/%d %s\n', i, numel(fileList), fileList{i});
        end
    end
end

R.files = fileList;
R.options = opts;
R.variants = opts.VariantNames;
R.sources = vertcat(res{:});
R.sources = vertcat(R.sources.sources);
R.events = vertcat(res{:});
R.events = vertcat(R.events.events);
R.summary = summarize(R.sources, R.events, opts.VariantNames);
R.bySize = bySize(R.events, opts.VariantNames, opts.SizeEdges);
R.seconds = toc(t0);
if opts.Verbose
    report(R, nV);
end
end

function opts = parseOptions(varargin)
ip = inputParser;
ip.addParameter('Variants', {struct()});
ip.addParameter('VariantNames', {});
ip.addParameter('TruthParams', struct(), @isstruct);
ip.addParameter('TruthScale', 1, @(x) isnumeric(x) && isscalar(x) && x >= 0);
ip.addParameter('InjectSizes', [], @(x) isnumeric(x) && all(x > 0));
ip.addParameter('InjectRate', 0.5, @(x) isnumeric(x) && isscalar(x) && x >= 0);
ip.addParameter('Seed', 0, @(x) isnumeric(x) && isscalar(x));
ip.addParameter('SizeEdges', [0 0.25 0.5 1 2 4 8 Inf], @(x) isnumeric(x) && issorted(x));
ip.addParameter('UseParallel', false, @(x) islogical(x) || isnumeric(x));
ip.addParameter('Verbose', true, @(x) islogical(x) || isnumeric(x));
ip.parse(varargin{:});
opts = ip.Results;
v = opts.Variants;
if isstruct(v)
    v = num2cell(v);
end
if ~iscell(v) || isempty(v) || ~all(cellfun(@isstruct, v))
    error('subproblemRefitTest:variants', 'Variants must be a struct, a struct array or a cell array of structs.');
end
opts.Variants = v(:)';
names = opts.VariantNames;
if isempty(names)
    names = arrayfun(@(k) sprintf('refit%d', k), 1:numel(v), 'UniformOutput', false);
end
names = cellstr(names);
if numel(names) ~= numel(v)
    error('subproblemRefitTest:variantNames', 'Give one name per variant.');
end
if numel(unique(names)) < numel(names) %the names group the results
    error('subproblemRefitTest:variantNames', 'Give each variant a name of its own.');
end
opts.VariantNames = names(:)';
opts.UseParallel = logical(opts.UseParallel);
opts.Verbose = logical(opts.Verbose);
end

function fileList = listFiles(files)
if (ischar(files) || (isstring(files) && isscalar(files))) && isfolder(files)
    d = dir(fullfile(char(files), '**', 'sub_*.mat'));
    fileList = fullfile({d.folder}, {d.name});
    fileList = sort(fileList(:));
elseif (ischar(files) || (isstring(files) && isscalar(files))) && endsWith(files, '.txt')
    fileList = strtrim(splitlines(string(fileread(files))));
    fileList = cellstr(fileList(fileList ~= ""));
else
    fileList = cellstr(files);
    fileList = fileList(:);
end
missing = fileList(~cellfun(@isfile, fileList));
if ~isempty(missing)
    error('subproblemRefitTest:noFiles', 'No such folder or subproblem file: %s', missing{1});
end
end

function p = withOverrides(p, overrides)
%the saved params with overrides; the fields extraction derives from them
%(tau_full, k, k2, the window lengths in samples, the footprint and surround
%filters) are recomputed when one they depend on is overridden, as SILo
%computes them, unless the override sets them too
fixed = intersect(fieldnames(overrides), {'photonScale', 'dXY', 'selRadius'});
if ~isempty(fixed)
    error('subproblemRefitTest:fixedParam', ...
        'params.%s is fixed by the saved subproblem (its data are normalized and its pixels selected with it).', fixed{1});
end
p = withFields(p, overrides);
base = {'tau_s', 'tau_full', 'tau2_s', 'analyzeHz', 'sigma_px', 'baselineWindow_Glu_s', 'denoiseWindow_s'};
if any(isfield(overrides, base))
    if ~isfield(overrides, 'tau_full') %an explicit tau_full sets the kernel
        p.tau_full = p.tau_s * p.analyzeHz;
    end
    p = setParamsExtractTrial(p);
    p = withFields(p, overrides);
end
end

function p = withFields(p, overrides)
for f = fieldnames(overrides)'
    p.(f{1}) = overrides.(f{1});
end
end

function out = oneSubproblem(file, i, opts)
d = load(file, 'Y_p', 'F_inv_p', 'sources_p', 'selPix_p', 'params');
Y = double(d.Y_p(:, :, 1));
F = double(d.F_inv_p);
p0 = withOverrides(d.params, opts.TruthParams);
[H0, S0, F00, ~, ~, ~, B0] = extractSubproblem(Y, F, d.sources_p, d.selPix_p, p0);
B0 = B0(:, :, 1); F00 = F00(:, :, 1);
S0 = opts.TruthScale * S0(:, :, 1);
X0 = convn(S0, p0.k, 'same');
nT = size(X0, 2); nS = size(H0, 2);
L = ceil(2.5 * p0.tau_full);
peakH = max(H0, [], 1);

prevRng = rng;
restoreRng = onCleanup(@() rng(prevRng));
rng(opts.Seed + i);
ev = cell(nS, 1);
for s = 1:nS
    ev{s} = findEvents(X0(s, :), L);
end
inj = cell(nS, 1);
if ~isempty(opts.InjectSizes) && opts.InjectRate > 0
    hz = p0.analyzeHz;
    nInj = round(opts.InjectRate * nT / hz);
    kPeak = max(p0.k);
    kc = floor(numel(p0.k)/2) + 1; %the kernel's centre: an event at t peaks at t
    kIx = kc + (-2:L); kIx = kIx(kIx >= 1 & kIx <= numel(p0.k)); %the window, within the kernel's support
    kWin = sum(p0.k(kIx)) / kPeak; %an event's mass in its window, per unit of peak
    for s = 1:nS
        [t, a] = injectionTimes(X0(s, :), ev{s}, nInj, L, opts.InjectSizes / peakH(s), kWin);
        S0(s, t) = S0(s, t) + a / kPeak;
        inj{s} = t;
    end
    X0 = convn(S0, p0.k, 'same');
end
T0 = max(B0 + H0 * X0, 0);
Z = F .* poissrnd(T0 ./ F); %mean T0, variance F.*T0: the fit's own noise model
clear restoreRng

nV = numel(opts.Variants);
srcRows = cell(nV, 1); evRows = cell(nV, 1);
for v = 1:nV
    pv = withOverrides(d.params, opts.Variants{v});
    [~, S1, F01, LS1] = extractSubproblem(Z, F, d.sources_p, d.selPix_p, pv);
    X1 = convn(S1(:, :, 1), pv.k, 'same');
    F01 = F01(:, :, 1); LS1 = LS1(:, :, 1);
    srcV = cell(nS, 1); evV = cell(nS, 1);
    for s = 1:nS
        x0 = X0(s, :); x1 = X1(s, :);
        act = x0 > 0.2 * max(x0);
        trueMass = sum(x0);
        if trueMass > 0
            ratios = [sum(x1) / trueMass, sum(x1 .* x0) / sum(x0.^2), sum(LS1(s, :)) / trueMass];
        else %no true activity
            ratios = nan(1, 3);
        end
        srcV{s} = table(string(file), s, string(opts.VariantNames{v}), trueMass, ratios(1), ratios(2), ...
            mean(F01(s, :)) / mean(F00(s, :)), mean(F01(s, act)) / mean(F00(s, act)), ratios(3), peakH(s), ...
            (H0(:, s)' * mean(B0, 2)) / sum(H0(:, s)), ...
            'VariableNames', {'file', 'source', 'variant', 'trueMass', 'recovered', 'slope', ...
            'F0ratio', 'F0ratioActive', 'lsRecovered', 'peakH', 'baseline'});
        e = [ev{s}(:); inj{s}(:)];
        isInj = [false(numel(ev{s}), 1); true(numel(inj{s}), 1)];
        m0 = zeros(numel(e), 1); m1 = zeros(numel(e), 1);
        for j = 1:numel(e)
            w = e(j) + (-2:L);
            m0(j) = sum(x0(w)); m1(j) = sum(x1(w));
        end
        evV{s} = table(repmat(string(file), numel(e), 1), repmat(s, numel(e), 1), ...
            repmat(string(opts.VariantNames{v}), numel(e), 1), reshape(x0(e), [], 1) * peakH(s), m0, m1, ...
            m1 > 0.25 * m0, isInj, 'VariableNames', ...
            {'file', 'source', 'variant', 'size', 'trueMass', 'refitMass', 'found', 'injected'});
    end
    srcRows{v} = vertcat(srcV{:});
    evRows{v} = vertcat(evV{:});
end
out.sources = vertcat(srcRows{:});
out.events = vertcat(evRows{:});
end

function pk = findEvents(x0, L)
%local maxima above 2% of the largest, at least L samples after the previous
%maximum, away from the ends
nT = numel(x0);
pk = find(x0(2:end-1) > 0.02 * max(x0) & x0(2:end-1) >= x0(1:end-2) & x0(2:end-1) > x0(3:end)) + 1;
pk = pk(pk > 5 & pk < nT - L - 3);
keep = true(size(pk));
for q = 2:numel(pk)
    if pk(q) - pk(q-1) < L
        keep(q) = false;
    end
end
pk = pk(keep);
end

function [t, a] = injectionTimes(x0, pk, n, L, peaks, kWin)
%up to n random times t for events of peak activity a (the values of peaks in
%turn): no true event within 2L+4 samples, less than 10% of the event's mass
%(a*kWin) already in its window, and 2L+4 samples apart
nT = numel(x0);
cand = 6:(nT - L - 4);
if ~isempty(pk)
    near = any(abs(cand(:) - pk(:).') <= 2*L + 4, 2);
    cand = cand(~near);
end
cand = cand(randperm(numel(cand)));
t = zeros(1, 0); a = zeros(1, 0);
for j = 1:numel(cand)
    if numel(t) >= n
        break
    end
    c = cand(j);
    aj = peaks(1 + mod(numel(t), numel(peaks)));
    if all(abs(c - t) > 2*L + 4) && sum(x0(c + (-2:L))) < 0.1 * aj * kWin
        t(end+1) = c; %#ok<AGROW>
        a(end+1) = aj; %#ok<AGROW>
    end
end
[t, order] = sort(t);
a = a(order);
end

function S = summarize(src, ev, names)
rows = cell(numel(names), 1);
for v = 1:numel(names)
    s = src(src.variant == names{v}, :);
    e = ev(ev.variant == names{v} & ~ev.injected, :);
    q = @(x, pct) prctile(x, pct); %prctile leaves out NaN
    rows{v} = table(string(names{v}), height(s), height(e), median(s.recovered, 'omitnan'), q(s.recovered, 25), q(s.recovered, 75), ...
        sum(e.refitMass) / sum(e.trueMass), median(s.slope, 'omitnan'), q(s.slope, 25), q(s.slope, 75), ...
        median(s.F0ratio, 'omitnan'), median(s.F0ratioActive, 'omitnan'), median(s.lsRecovered, 'omitnan'), ...
        'VariableNames', {'variant', 'nSources', 'nEvents', 'recovered', 'recoveredQ1', 'recoveredQ3', ...
        'pooled', 'slope', 'slopeQ1', 'slopeQ3', 'F0ratio', 'F0ratioActive', 'lsRecovered'});
end
S = vertcat(rows{:});
end

function B = bySize(ev, names, edges)
rows = {};
for v = 1:numel(names)
    for injected = [false true]
        e = ev(ev.variant == names{v} & ev.injected == injected, :);
        for b = 1:numel(edges) - 1
            m = e.size >= edges(b) & e.size < edges(b+1);
            if ~any(m)
                continue
            end
            f = m & e.found; %recoveredFound is NaN when none is found (0/0)
            rows{end+1, 1} = table(string(names{v}), injected, edges(b), edges(b+1), nnz(m), nnz(f) / nnz(m), ...
                sum(e.refitMass(m)) / sum(e.trueMass(m)), sum(e.refitMass(f)) / sum(e.trueMass(f)), ...
                'VariableNames', {'variant', 'injected', 'sizeFrom', 'sizeTo', 'n', 'shareFound', ...
                'recoveredAll', 'recoveredFound'}); %#ok<AGROW>
        end
    end
end
if isempty(rows) %no events: an empty table with the same columns
    B = table('Size', [0 8], 'VariableTypes', {'string', 'logical', 'double', 'double', 'double', 'double', 'double', 'double'}, ...
        'VariableNames', {'variant', 'injected', 'sizeFrom', 'sizeTo', 'n', 'shareFound', 'recoveredAll', 'recoveredFound'});
else
    B = vertcat(rows{:});
end
end

function report(R, nV)
o = R.options;
fprintf('\n%d subproblems, %.0f s. Truth scale %.2f', numel(R.files), R.seconds, o.TruthScale);
if ~isempty(o.InjectSizes)
    fprintf(', injected sizes %s at %.2f per source and second', mat2str(o.InjectSizes), o.InjectRate);
end
fprintf('\n');
for v = 1:nV
    s = R.summary(v, :);
    src = R.sources(R.sources.variant == s.variant, :);
    fprintf('\n== %s (overrides: %s)\n', s.variant, strjoin(fieldnames(o.Variants{v})', ', '));
    fprintf('%d sources, %d true events. Baseline under the footprints: median %.2f (quartiles %.2f-%.2f) per pixel and sample\n', ...
        s.nSources, s.nEvents, median(src.baseline), prctile(src.baseline, 25), prctile(src.baseline, 75));
    fprintf('Activity mass recovered per source (refit / truth): median %.2f, quartiles %.2f-%.2f; pooled %.2f\n', ...
        s.recovered, s.recoveredQ1, s.recoveredQ3, s.pooled);
    fprintf('Slope of the refit on the true activity (through 0): median %.2f, quartiles %.2f-%.2f\n', s.slope, s.slopeQ1, s.slopeQ3);
    fprintf('F0 refit / true: all samples median %.3f; during activity median %.3f\n', s.F0ratio, s.F0ratioActive);
    fprintf('Least-squares trace mass / true activity mass: median %.2f\n', s.lsRecovered);
    b = R.bySize(R.bySize.variant == s.variant, :);
    for injected = [false true]
        bi = b(b.injected == injected, :);
        if isempty(bi), continue, end
        if injected
            fprintf('Injected events:\n');
        else
            fprintf('By true event size (peak-pixel amplitude, normalized units):\n');
        end
        fprintf('   size          events   share found   recovered (all)   recovered (found)\n');
        for k = 1:height(bi)
            fprintf('   %4.2f-%-6.2f %8d   %8.2f      %8.2f          %8.2f\n', bi.sizeFrom(k), bi.sizeTo(k), ...
                bi.n(k), bi.shareFound(k), bi.recoveredAll(k), bi.recoveredFound(k));
        end
    end
end
end
