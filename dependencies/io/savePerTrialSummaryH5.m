function savePerTrialSummaryH5(filename, exptSummary, trialTable)
%SAVEPERTRIALSUMMARYH5 Write SILo per-trial summary to HDF5.
%
%   savePerTrialSummaryH5(FILENAME, EXPTSUMMARY, TRIALTABLE) writes the
%   per-trial summary produced by SILo.m to disk using the per-path tree
%   documented in the README.md "Experiment Summary" section:
%       /Path{n}/sources/temporal/per_trial_SNR
%       /Path{n}/sources/spatial/{per_trial_profiles, per_trial_coords}
%       /Path{n}/visualizations/{per_trial_mean_im, per_trial_F0_im,
%                                per_trial_act_im, per_trial_act_im_peaks,
%                                per_trial_num_peaks}
%       /Path{n}/user_rois/{per_trial_F0, per_trial_F0_coverage}
%
%   The trial axis matches trial_table.h5 (all analysis trials). Trials
%   without source extraction or alignment data are left as NaN in the
%   corresponding arrays.

if nargin < 3 || isempty(filename) || isempty(exptSummary) || isempty(trialTable)
    error('savePerTrialSummaryH5:MissingInput', ...
        'filename, exptSummary, and trialTable are required.');
end

nDMDs = size(trialTable.filename, 1);
nTrials = size(trialTable.filename, 2);

if ~isfield(exptSummary, 'E') || isempty(exptSummary.E)
    exptSummary.E = cell(nTrials, nDMDs);
elseif size(exptSummary.E, 1) < nTrials || size(exptSummary.E, 2) < nDMDs
    exptSummary.E{nTrials, nDMDs} = [];
end

s = struct();
for DMDix = 1:nDMDs
    if isempty(exptSummary.meanIM{DMDix})
        continue
    end
    pathName = sprintf('Path%d', DMDix);
    s.(pathName) = buildPerTrialPathGroup(exptSummary, DMDix, nTrials);
end

saveStructToH5(s, filename);
end


function ps = buildPerTrialPathGroup(exptSummary, DMDix, nTrials)
ps = struct();
selPix = exptSummary.selPix{DMDix};
[imRows, imCols] = size(selPix);

ps.visualizations = buildPerTrialVisualizations(exptSummary, DMDix, nTrials);

% per-trial background (F0) images: rows x cols x channels x trials, original
% channel order; empty if no trial has source extraction results
f0Stack = buildPerTrialF0Stack(exptSummary, DMDix, nTrials, selPix);
if ~isempty(f0Stack)
    ps.visualizations.per_trial_F0_im = permute(reshape(f0Stack, ...
        imRows, imCols, size(f0Stack, 3), 1, nTrials), [5 3 4 1 2]);
    ur = buildPerTrialUserRois(exptSummary, DMDix, f0Stack);
    if ~isempty(fieldnames(ur))
        ps.user_rois = ur;
    end
end
src = buildPerTrialSources(exptSummary, DMDix, nTrials, selPix, imRows, imCols);
if ~isempty(fieldnames(src))
    ps.sources = src;
end
end


function vis = buildPerTrialVisualizations(exptSummary, DMDix, nTrials)
vis = struct();

meanAligned = exptSummary.perTrialMeanIMsAligned{DMDix};
actAligned = exptSummary.perTrialActIMsAligned{DMDix};

H = size(meanAligned, 1);
W = size(meanAligned, 2);
C = size(meanAligned, 3);

% Always write the full trial_table trial axis. Invalid / skipped trials
% (all-NaN actAligned from SILo) are left as all-NaN slices. meanAligned
% is not preallocated in SILo, so skipped slots can be zero rather than NaN.
perTrialMeanIm = nan(nTrials, C, 1, H, W, 'like', meanAligned);
for tIx = 1:nTrials
    if all(isnan(actAligned(:, :, 1, tIx)), 'all')
        continue
    end
    perTrialMeanIm(tIx, :, 1, :, :) = permute(meanAligned(:, :, :, tIx), [3 1 2]);
end
vis.per_trial_mean_im = perTrialMeanIm;

% actAligned is preallocated to nTrials with NaN for skipped trials in SILo.
vis.per_trial_act_im = permute(actAligned, [4 3 1 2]);

peaksCell = exptSummary.peaks{DMDix};
if ~isempty(peaksCell)
    actImPeaks = buildPerTrialActImPeaks(peaksCell, nTrials, exptSummary, DMDix);
    vis.per_trial_num_peaks = actImPeaks.per_trial_num_peaks;
    if isfield(actImPeaks, 'per_trial_act_im_peaks')
        vis.per_trial_act_im_peaks = actImPeaks.per_trial_act_im_peaks;
    end
end
end


function peaksOut = buildPerTrialActImPeaks(peaksCell, nTrials, exptSummary, DMDix)
% trials x max_peaks x 3 [z_loc, y_loc, x_loc], NaN-padded; z_loc is 0.
peaksOut = struct();

numPeaks = zeros(nTrials, 1);
for tIx = 1:nTrials
    if ~isfield(peaksCell{tIx}, 'row') || isempty(peaksCell{tIx}.row)
        continue
    end
    numPeaks(tIx) = numel(peaksCell{tIx}.row(:));
end
peaksOut.per_trial_num_peaks = int32(numPeaks);

maxPeaks = max(numPeaks);
if maxPeaks == 0
    return
end

motOutput = [];
if isfield(exptSummary, 'perTrialAlignmentOffsets') ...
        && numel(exptSummary.perTrialAlignmentOffsets) >= DMDix
    motOutput = exptSummary.perTrialAlignmentOffsets{DMDix};
end
warp = []; % common frame -> trial image, see trialWarpCoords
if isfield(exptSummary, 'perTrialAlignmentWarp') ...
        && numel(exptSummary.perTrialAlignmentWarp) >= DMDix
    warp = exptSummary.perTrialAlignmentWarp{DMDix};
elseif isfield(exptSummary, 'perTrialAlignmentTforms') ...
        && numel(exptSummary.perTrialAlignmentTforms) >= DMDix ...
        && ~isempty(exptSummary.perTrialAlignmentTforms{DMDix})
    warp = struct('affine', exptSummary.perTrialAlignmentTforms{DMDix}, 'nodes', zeros(0, 2)); % older summaries: affine maps only
end

perTrialActImPeaks = nan(nTrials, maxPeaks, 3);
for tIx = 1:nTrials
    nP = numPeaks(tIx);
    if nP == 0
        continue
    end
    pRows = double(peaksCell{tIx}.row(:));
    pCols = double(peaksCell{tIx}.col(:));
    %peaks are found in the trial's own registered image; map them into the common frame
    if ~isempty(warp) && size(warp.affine, 3) >= tIx && all(~isnan(warp.affine(:, :, tIx)), 'all')
        [pRows, pCols] = trialWarpInverse(warp, tIx, pRows, pCols);
    elseif ~isempty(motOutput) && size(motOutput, 2) >= tIx
        pRows = pRows - motOutput(1, tIx);
        pCols = pCols - motOutput(2, tIx);
    end
    perTrialActImPeaks(tIx, 1:nP, 1) = 0;
    perTrialActImPeaks(tIx, 1:nP, 2) = pRows - 1;
    perTrialActImPeaks(tIx, 1:nP, 3) = pCols - 1;
end
peaksOut.per_trial_act_im_peaks = perTrialActImPeaks;
end


function f0Stack = buildPerTrialF0Stack(exptSummary, DMDix, nTrials, selPix)
% rows x cols x channels x trials; trials without F0 are all-NaN.
f0Stack = [];
params = exptSummary.params;
numChannels = params.numChannels;
activityChannel = params.activityChannel;
orderedChannels = [activityChannel:numChannels, 1:activityChannel-1];
for tIx = 1:nTrials
    img = buildF0Image(exptSummary.E{tIx, DMDix}, selPix, numChannels, orderedChannels);
    if isempty(img)
        continue
    end
    if isempty(f0Stack)
        f0Stack = nan(size(img, 1), size(img, 2), numChannels, nTrials, 'single');
    end
    f0Stack(:, :, :, tIx) = img;
end
end


function ur = buildPerTrialUserRois(exptSummary, DMDix, f0Stack)
% Mean F0 over each user ROI mask per trial (rois x channels x trials), and
% the fraction of mask pixels with an F0 estimate in the activity channel
% (rois x trials). Masks are in the same across-trial aligned frame as F0.
ur = struct();
roiData = exptSummary.userROIs{DMDix};
nRois = numel(roiData);
if nRois == 0
    return
end
[imRows, imCols, numChannels, nTrials] = size(f0Stack);
activityChannel = exptSummary.params.activityChannel;
f0Px = reshape(f0Stack, imRows*imCols, numChannels, nTrials);

perTrialF0 = nan(nRois, numChannels, nTrials);
coverage = zeros(nRois, nTrials);
for rIx = 1:nRois
    if ~isfield(roiData{rIx}, 'mask') || isempty(roiData{rIx}.mask)
        continue
    end
    tmp = logical(roiData{rIx}.mask);
    mask = false(imRows, imCols);
    nr = min(imRows, size(tmp, 1)); nc = min(imCols, size(tmp, 2));
    mask(1:nr, 1:nc) = tmp(1:nr, 1:nc);
    if ~any(mask(:))
        continue
    end
    vals = f0Px(mask(:), :, :); % mask pixels x channels x trials
    perTrialF0(rIx, :, :) = mean(vals, 1, 'omitnan');
    coverage(rIx, :) = reshape(mean(~isnan(vals(:, activityChannel, :)), 1), 1, nTrials);
end
ur.per_trial_F0 = perTrialF0;
ur.per_trial_F0_coverage = coverage;
end


function src = buildPerTrialSources(exptSummary, DMDix, nTrials, selPix, imRows, imCols)
src = struct();

nSources = 0;
for tIx = 1:nTrials
    E = exptSummary.E{tIx, DMDix};
    if ~isempty(E) && isfield(E, 'footprints') && ~isempty(E.footprints)
        nSources = size(E.footprints, 2);
        break
    end
end
if nSources == 0
    return
end

selFlat = find(selPix);
[rPix, cPix] = ind2sub([imRows, imCols], selFlat);

perTrialProfiles = nan(nTrials, nSources, 1, imRows, imCols, 'single');
perTrialCoords = nan(nTrials, nSources, 3);
perTrialCoords(:, :, 1) = 0;
perTrialSnr = nan(nTrials, nSources);

for tIx = 1:nTrials
    E = exptSummary.E{tIx, DMDix};
    if isempty(E) || ~isfield(E, 'footprints') || isempty(E.footprints)
        continue
    end
    H = E.footprints;
    if size(H, 2) ~= nSources
        continue
    end

    perTrialSnr(tIx, :) = double(E.SNR(:));

    for sIx = 1:nSources
        w = double(H(:, sIx));
        w(isnan(w)) = 0;
        img = zeros(imRows, imCols, 'single');
        img(selPix) = w;
        perTrialProfiles(tIx, sIx, 1, :, :) = img;
        wsum = sum(w, 'omitnan');
        if wsum > 0
            perTrialCoords(tIx, sIx, 2) = sum(w .* rPix, 'omitnan') / wsum - 1;
            perTrialCoords(tIx, sIx, 3) = sum(w .* cPix, 'omitnan') / wsum - 1;
        end
    end
end

sp = struct();
sp.per_trial_profiles = perTrialProfiles;
sp.per_trial_coords = perTrialCoords;
src.spatial = sp;

temp = struct();
temp.per_trial_SNR = perTrialSnr;
src.temporal = temp;
end
