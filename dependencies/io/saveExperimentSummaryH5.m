function saveExperimentSummaryH5(filename, exptSummary, trialTable)
%SAVEEXPERIMENTSUMMARYH5 Write SILo experiment summary to HDF5.
%
%   saveExperimentSummaryH5(FILENAME, EXPTSUMMARY, TRIALTABLE) writes the
%   experiment summary produced by SILo.m to disk using the per-path tree
%   documented in the README.md "Experiment Summary" section:
%       /params
%       /Path{n}/Z_depths
%       /Path{n}/frame_info/{offline,online}{X,Y,Z}shifts,
%                          trial_num_frames, frame_line_idxs, discard_frames
%       /Path{n}/visualizations/{mean_im, act_im, act_im_peaks}
%       /Path{n}/global/F
%       /Path{n}/user_rois/{labels, mask, F, Fsvd}
%       /Path{n}/sources/spatial/{profiles, coords}
%       /Path{n}/sources/temporal/{dF_ls, dF_denoised, events, F0, SNR}
%
%   Per-trial traces are concatenated along the time axis for the valid
%   trials in exptSummary.E; invalid (skipped) trials contribute 0 frames
%   to the "total frames" axis. The SLAP2 reference stack is not duplicated
%   here; it remains in trial_table.h5 under slap2_info/ref_stack.

if nargin < 3 || isempty(filename) || isempty(exptSummary) || isempty(trialTable)
    error('saveExperimentSummaryH5:MissingInput', ...
        'filename, exptSummary, and trialTable are required.');
end

params = exptSummary.params;
isSlap2 = logical(params.isSLAP2);
nDMDs = size(trialTable.filename, 1);
nTrials = size(trialTable.filename, 2);

% Channel layout. Internally the pipeline stores per-channel data in
% "ordered" order (activity channel first, then the rest cycling around).
% The H5 readout uses original channel order so downstream consumers can
% rely on params.activityChannel to find the glutamate channel.
numChannels = params.numChannels;
activityChannel = params.activityChannel;
orderedChannels = [activityChannel:numChannels, 1:activityChannel-1];

% Pad exptSummary.E to (nTrials, nDMDs) so direct cell indexing is safe.
% SILo.m only assigns exptSummary.E(:,DMDix) when sources were found, so it
% can be missing entirely or sized smaller than (nTrials, nDMDs).
if ~isfield(exptSummary, 'E') || isempty(exptSummary.E)
    exptSummary.E = cell(nTrials, nDMDs);
elseif size(exptSummary.E, 1) < nTrials || size(exptSummary.E, 2) < nDMDs
    exptSummary.E{nTrials, nDMDs} = [];
end

s = struct();
s.params = params;

for DMDix = 1:nDMDs
    if isempty(exptSummary.meanIM{DMDix})
        % Every trial for this DMD failed alignment; nothing to save.
        continue
    end
    pathName = sprintf('Path%d', DMDix);
    s.(pathName) = buildPathGroup(exptSummary, DMDix, isSlap2, nTrials, ...
        numChannels, orderedChannels);
end

saveStructToH5(s, filename);
end


function ps = buildPathGroup(exptSummary, DMDix, isSlap2, nTrials, ...
    numChannels, orderedChannels)
ps = struct();

[framesPerTrial, totalFrames] = collectFrameCounts(exptSummary, DMDix, nTrials);
selPix = exptSummary.selPix{DMDix};
[imRows, imCols] = size(selPix);

% Z_depths (fastz x 1) — SLAP2 only per README, but written if available.
zVal = double(exptSummary.Z(DMDix));
if isSlap2 && ~isnan(zVal)
    ps.Z_depths = zVal(:);
end

ps.frame_info = buildFrameInfo(exptSummary, DMDix, isSlap2, ...
    framesPerTrial, totalFrames, nTrials);
ps.visualizations = buildVisualizations(exptSummary, DMDix);
ps.global = buildGlobal(exptSummary, DMDix, framesPerTrial, totalFrames, nTrials);
ps.user_rois = buildUserRois(exptSummary, DMDix, framesPerTrial, totalFrames, ...
    nTrials, imRows, imCols, numChannels, orderedChannels);
ps.sources = buildSources(exptSummary, DMDix, framesPerTrial, totalFrames, ...
    nTrials, selPix, numChannels, orderedChannels);
end


function [framesPerTrial, totalFrames] = collectFrameCounts(exptSummary, DMDix, nTrials)
framesPerTrial = zeros(nTrials, 1);
for tIx = 1:nTrials
    E = exptSummary.E{tIx, DMDix};
    if isempty(E)
        continue
    end
    framesPerTrial(tIx) = trialFrameCount(E);
end
totalFrames = sum(framesPerTrial);
end


function n = trialFrameCount(E)
n = 0;
if isfield(E, 'dF') && isfield(E.dF, 'ls') && ~isempty(E.dF.ls)
    n = size(E.dF.ls, 2);
elseif isfield(E, 'discardFrames')
    n = numel(E.discardFrames);
elseif isfield(E, 'global') && isfield(E.global, 'F') && ~isempty(E.global.F)
    n = size(E.global.F, 2);
end
end


function fi = buildFrameInfo(exptSummary, DMDix, isSlap2, framesPerTrial, totalFrames, nTrials)
fi = struct();
fi.trial_num_frames = int32(framesPerTrial(:));

offlineX = nan(totalFrames, 1);
offlineY = nan(totalFrames, 1);
offlineZ = nan(totalFrames, 1);
onlineX = nan(totalFrames, 1);
onlineY = nan(totalFrames, 1);
onlineZ = nan(totalFrames, 1);
frameLineIdxs = zeros(totalFrames, 1);
discard = false(totalFrames, 1);

motOffsets = [];
if isfield(exptSummary, 'perTrialAlignmentOffsets') && numel(exptSummary.perTrialAlignmentOffsets) >= DMDix
    motOffsets = exptSummary.perTrialAlignmentOffsets{DMDix};
end

frameOff = 0;
for tIx = 1:nTrials
    n = framesPerTrial(tIx);
    if n == 0
        continue
    end
    idxs = frameOff + (1:n);
    E = exptSummary.E{tIx, DMDix};

    if isfield(E, 'discardFrames')
        discard(idxs) = logical(E.discardFrames(:));
    end

    aData = exptSummary.aData{tIx, DMDix};
    dr = 0; dc = 0;
    if ~isempty(motOffsets) && size(motOffsets, 2) >= tIx
        dr = motOffsets(1, tIx);
        dc = motOffsets(2, tIx);
    end

    frameLines = double(E.frameLines(:))';
    frameLineIdxs(idxs) = frameLines;

    if isSlap2
        if isfield(aData, 'DSframes') && ~isempty(aData.DSframes)
            DSf = double(aData.DSframes(:))';
            if isfield(aData, 'motionDSc')
                offlineX(idxs) = upsampleMotion(DSf, double(aData.motionDSc(:))', frameLines, 'pchip') + dc;
            end
            if isfield(aData, 'motionDSr')
                offlineY(idxs) = upsampleMotion(DSf, double(aData.motionDSr(:))', frameLines, 'pchip') + dr;
            end
            if isfield(aData, 'motionDSz')
                offlineZ(idxs) = upsampleMotion(DSf, double(aData.motionDSz(:))', frameLines, 'pchip');
            end
            slap2 = struct();
            if isfield(aData, 'slap2'), slap2 = aData.slap2; end
            if isfield(slap2, 'onlineMotionXshift')
                onlineX(idxs) = upsampleMotion(DSf, double(slap2.onlineMotionXshift(:))', frameLines, 'nearest');
            end
            if isfield(slap2, 'onlineMotionYshift')
                onlineY(idxs) = upsampleMotion(DSf, double(slap2.onlineMotionYshift(:))', frameLines, 'nearest');
            end
            if isfield(slap2, 'onlineMotionZshift')
                onlineZ(idxs) = upsampleMotion(DSf, double(slap2.onlineMotionZshift(:))', frameLines, 'nearest');
            end
        end
    else
        if isfield(aData, 'motionC') && ~isempty(aData.motionC)
            offlineX(idxs) = double(aData.motionC(:)) + dc;
        end
        if isfield(aData, 'motionR') && ~isempty(aData.motionR)
            offlineY(idxs) = double(aData.motionR(:)) + dr;
        end
        if isfield(aData, 'motionZ') && ~isempty(aData.motionZ)
            offlineZ(idxs) = double(aData.motionZ(:));
        end
    end

    frameOff = frameOff + n;
end

fi.offlineXshifts = offlineX;
fi.offlineYshifts = offlineY;
if ~all(isnan(offlineZ)), fi.offlineZshifts = offlineZ; end
fi.frame_line_idxs = int32(frameLineIdxs);
fi.discard_frames = discard;
if isSlap2
    fi.onlineXshifts = onlineX;
    fi.onlineYshifts = onlineY;
    fi.onlineZshifts = onlineZ;
end
end


function v = upsampleMotion(DSf, motionDS, frameLines, method)
%motion at each frame line from the downsampled motion; constant when a
%(short) trial has a single downsampled frame, which interp1 cannot use
if numel(DSf) > 1
    v = interp1(DSf, motionDS, frameLines, method, 'extrap');
else
    v = repmat(motionDS(1), size(frameLines));
end
end


function vis = buildVisualizations(exptSummary, DMDix)
vis = struct();

% rows x cols x channels -> channels x fastz(1) x rows x cols
mIm = exptSummary.meanIM{DMDix};
H = size(mIm, 1); W = size(mIm, 2); C = size(mIm, 3);
mIm3 = reshape(mIm, H, W, C);
vis.mean_im = reshape(permute(mIm3, [3 1 2]), C, 1, H, W);

% rows x cols x 1 -> fastz(1) x rows x cols
aIm = exptSummary.actIM{DMDix};
vis.act_im = reshape(aIm(:, :, 1), 1, size(aIm, 1), size(aIm, 2));

% sources x 3 [z_loc, y_loc, x_loc], 0-indexed; z_loc fixed at 0.
sources = exptSummary.sources{DMDix};
if isfield(sources, 'R') && isfield(sources, 'C') && ~isempty(sources.R)
    nSources = numel(sources.R);
    actImPeaks = zeros(nSources, 3);
    actImPeaks(:, 2) = double(sources.R(:)) - 1;
    actImPeaks(:, 3) = double(sources.C(:)) - 1;
    vis.act_im_peaks = actImPeaks;
end
end


function g = buildGlobal(exptSummary, DMDix, framesPerTrial, totalFrames, nTrials)
g = struct();
if totalFrames == 0
    return
end
nChannels = 0;
for tIx = 1:nTrials
    E = exptSummary.E{tIx, DMDix};
    if ~isempty(E) && isfield(E, 'global') && isfield(E.global, 'F') && ~isempty(E.global.F)
        nChannels = size(E.global.F, 1);
        break
    end
end
if nChannels == 0
    return
end
F = nan(nChannels, totalFrames);
frameOff = 0;
for tIx = 1:nTrials
    n = framesPerTrial(tIx);
    if n == 0, continue, end
    E = exptSummary.E{tIx, DMDix};
    if isfield(E, 'global') && isfield(E.global, 'F') && ~isempty(E.global.F)
        F(:, frameOff + (1:n)) = E.global.F;
    end
    frameOff = frameOff + n;
end
g.F = F;
end


function ur = buildUserRois(exptSummary, DMDix, framesPerTrial, totalFrames, ...
    nTrials, imRows, imCols, numChannels, orderedChannels)
ur = struct();
roiData = exptSummary.userROIs{DMDix};
nRois = numel(roiData);
if nRois == 0
    return
end

labels = strings(nRois, 1);
masks = false(nRois, 1, imRows, imCols);
for rIx = 1:nRois
    if isfield(roiData{rIx}, 'Label')
        labels(rIx) = string(roiData{rIx}.Label);
    end
    if isfield(roiData{rIx}, 'mask') && ~isempty(roiData{rIx}.mask)
        masks(rIx, 1, :, :) = logical(roiData{rIx}.mask);
    end
end
ur.labels = cellstr(labels);
ur.mask = uint8(masks);

if totalFrames == 0
    return
end
% E.ROIs.F/Fsvd are (nRois x numChannels x frames) in ordered channel order;
% remap to original channel order on insertion.
F = nan(nRois, numChannels, totalFrames);
Fsvd = nan(nRois, numChannels, totalFrames);
frameOff = 0;
for tIx = 1:nTrials
    n = framesPerTrial(tIx);
    if n == 0, continue, end
    E = exptSummary.E{tIx, DMDix};
    if isfield(E, 'ROIs')
        if isfield(E.ROIs, 'F') && ~isempty(E.ROIs.F)
            F(:, orderedChannels, frameOff + (1:n)) = E.ROIs.F;
        end
        if isfield(E.ROIs, 'Fsvd') && ~isempty(E.ROIs.Fsvd)
            Fsvd(:, orderedChannels, frameOff + (1:n)) = E.ROIs.Fsvd;
        end
    end
    frameOff = frameOff + n;
end
ur.F = F;
ur.Fsvd = Fsvd;
end


function src = buildSources(exptSummary, DMDix, framesPerTrial, totalFrames, ...
    nTrials, selPix, numChannels, orderedChannels)
src = struct();

% Find first valid trial to grab spatial footprints and shape info
firstValid = 0;
E0 = [];
for tIx = 1:nTrials
    Ec = exptSummary.E{tIx, DMDix};
    if ~isempty(Ec) && isfield(Ec, 'footprints') && ~isempty(Ec.footprints)
        firstValid = tIx;
        E0 = Ec;
        break
    end
end
if firstValid == 0
    return
end
nSources = size(E0.footprints, 2);
% extractSources analyses at most 2 channels ([activity, secondary]); the
% size of the channel axis of dF therefore reports the analyzed count, not
% the recording's total channel count.
nAnalyzed = 0;
if isfield(E0, 'dF') && isfield(E0.dF, 'ls') && ~isempty(E0.dF.ls)
    nAnalyzed = size(E0.dF.ls, 3);
end

% spatial: average footprint images across trials, then centroid of each
% averaged profile (y_loc = row, x_loc = column).
sp = struct();
if ~isempty(selPix)
    [imRows, imCols] = size(selPix);
    profileSum = zeros(nSources, 1, imRows, imCols);
    nProfileTrials = 0;
    for tIx = 1:nTrials
        E = exptSummary.E{tIx, DMDix};
        if isempty(E) || ~isfield(E, 'footprints') || isempty(E.footprints)
            continue
        end
        nProfileTrials = nProfileTrials + 1;
        for sIx = 1:nSources
            w = double(E.footprints(:, sIx));
            w(isnan(w)) = 0;
            img = zeros(imRows, imCols, 'single');
            img(selPix) = w;
            profileSum(sIx, 1, :, :) = profileSum(sIx, 1, :, :) + reshape(img, 1, 1, imRows, imCols);
        end
    end
    profiles = single(profileSum / nProfileTrials);

    coords = zeros(nSources, 3);
    % coords are 0-indexed for HDF5/Python interop: [z_loc, y_loc, x_loc]
    % matching image axis order (fastz, rows, cols). z_loc is the 0-based
    % index into the fastz axis of profiles (always 0 for currently
    % supported single-plane recordings); y_loc/x_loc are row/column
    % centroids in the [0, dim-1] convention.
    [rGrid, cGrid] = ndgrid(1:imRows, 1:imCols);
    for sIx = 1:nSources
        img = squeeze(profiles(sIx, 1, :, :));
        wsum = sum(img(:));
        if wsum > 0
            coords(sIx, 2) = sum(img(:) .* rGrid(:)) / wsum - 1;
            coords(sIx, 3) = sum(img(:) .* cGrid(:)) / wsum - 1;
        end
    end
    sp.profiles = profiles;
    sp.coords = coords;
end
src.spatial = sp;

% temporal: concatenate per-trial traces along total-frames axis. Traces
% come out of extractSources in ordered channel order with only the
% analyzed channel(s); write into the full numChannels axis (original
% order) so un-analyzed channels are left as NaN.
if totalFrames == 0 || nAnalyzed == 0
    return
end
analyzedOriginal = orderedChannels(1:nAnalyzed);
dF_ls = nan(nSources, numChannels, totalFrames);
dF_denoised = nan(nSources, numChannels, totalFrames);
events = nan(nSources, numChannels, totalFrames);
F0 = nan(nSources, numChannels, totalFrames);

frameOff = 0;
for tIx = 1:nTrials
    n = framesPerTrial(tIx);
    if n == 0, continue, end
    E = exptSummary.E{tIx, DMDix};
    idxs = frameOff + (1:n);
    if isfield(E, 'dF')
        if isfield(E.dF, 'ls') && ~isempty(E.dF.ls)
            dF_ls(:, analyzedOriginal, idxs) = permute(E.dF.ls, [1 3 2]);
        end
        if isfield(E.dF, 'denoised') && ~isempty(E.dF.denoised)
            dF_denoised(:, analyzedOriginal, idxs) = permute(E.dF.denoised, [1 3 2]);
        end
        if isfield(E.dF, 'events') && ~isempty(E.dF.events)
            events(:, analyzedOriginal, idxs) = permute(E.dF.events, [1 3 2]);
        end
    end
    if isfield(E, 'F0') && ~isempty(E.F0)
        F0(:, analyzedOriginal, idxs) = permute(E.F0, [1 3 2]);
    end
    frameOff = frameOff + n;
end

temp = struct();
temp.dF_ls = dF_ls;
temp.dF_denoised = dF_denoised;
temp.events = events;
temp.F0 = F0;
snrAcrossTrials = nan(nSources, nTrials);
for tIx = 1:nTrials
    E = exptSummary.E{tIx, DMDix};
    if ~isempty(E) && isfield(E, 'SNR') && numel(E.SNR) == nSources
        snrAcrossTrials(:, tIx) = double(E.SNR(:));
    end
end
if any(~isnan(snrAcrossTrials), 'all')
    temp.SNR = mean(snrAcrossTrials, 2, 'omitnan');
end
src.temporal = temp;
end

