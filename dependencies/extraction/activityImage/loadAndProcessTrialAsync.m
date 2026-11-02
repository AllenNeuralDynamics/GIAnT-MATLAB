function  [IM, meanIM, IMc, aData, peaks, discardFrames, actAcc, halfIMs]= loadAndProcessTrialAsync(mocoDr, fn, numChannels, params)
% mocoDr: directory with motion-corrected outputs (downsampled registered TIFF/H5 and *_ALIGNMENTDATA.h5).
% halfIMs: rows x cols x channels x 2, the mean image of two independent
% halves of the frames (alternating 0.5 s blocks), for the reliability of the
% across-trial alignment (alignTrialsPatchwise).
    if endsWith(fn, '.h5')
        desc = h5info([mocoDr filesep fn]);
        IM = h5read([mocoDr filesep fn], ['/', desc.Datasets.Name]);
    else
        IM = ScanImageTiffWrapper([mocoDr filesep fn]);
    end
    IM = reshape(IM, size(IM,1), size(IM,2), numChannels, []); %deinterleave;

    meanIM= mean(IM,4, 'omitnan');
    nanPx = repmat(mean(isnan(IM), [3 4])>params.nanThresh, 1,1,size(IM,3));
    meanIM(nanPx) = nan;

    %halves from alternating blocks rather than alternating frames, because
    %neighbouring downsampled frames share noise; summed frame by frame to
    %avoid copying the movie
    blockLen = max(1, round(0.5*params.alignHz));
    half = mod(floor((0:size(IM,4)-1)/blockLen), 2) + 1;
    halfIMs = nan([size(meanIM, [1 2 3]) 2], 'like', meanIM);
    for h = 1:2
        frames = find(half == h);
        s = zeros(size(meanIM, [1 2 3])); n = zeros(size(meanIM, [1 2 3]));
        for f = frames
            x = double(IM(:,:,:,f));
            ok = ~isnan(x);
            x(~ok) = 0;
            s = s + x;
            n = n + ok;
        end
        m = s ./ n;
        m(repmat(1 - mean(n, 3)/numel(frames) > params.nanThresh, 1, 1, size(IM,3))) = nan;
        halfIMs(:,:,:,h) = m;
    end

    IM = squeeze(IM(:,:,params.activityChannel,:));

    %load alignment data
    fnStemEnd = strfind(fn, '_REGISTERED') -1;
    aData = loadStructFromH5([mocoDr filesep fn(1:fnStemEnd) '_ALIGNMENTDATA.h5']);
    aData.dsFac = 1; %SLAP2 data does not have downsampling per se
    params.dsFac = aData.dsFac;
    params.alignHz = aData.alignHz;
    nInitFrames = ceil(params.discardInitial_s * params.alignHz); %initial frames to discard

    %discard motion frames
    tmp = aData.recNegErr(:)- medfilt1(aData.recNegErr(:), round(4*params.alignHz)); %smoothExp(aData.recNegErr(:),'movmedian', ceil(2/(aData.frametime*aData.dsFac))); %-smoothdata(aData.aRankCorrDS,2, 'movmedian', ceil(2/aData.frametime));
    tmp = -tmp./(min(-0.005, prctile(tmp,5))); %normalize to the median-to-5th prctile interval, or 0.5% of intensity, whichever is larger
    thresh = params.motionThresh; %decrease thresh to be more stringent on motion correction
    window = 2*ceil(0.025*params.alignHz)+1;% a window in time to censor aronud motion events, ~50 milliseconds;
    discardFrames = imclose(imdilate(tmp>thresh, ones(window,1)) & (tmp>(thresh/2)), ones(window,1));
    discardFrames(1:nInitFrames) = true;
    IM(:,:,discardFrames) = nan;
    if isfield(aData, 'slap2') && isfield(aData.slap2, 'varFacDS')
        vIM = aData.slap2.varFacDS;
        aData.slap2 = rmfield(aData.slap2, 'varFacDS');
    else
        vIM = [];
    end

    try
       [IMc, peaks, actAcc] = localizeSources_vIM(IM, vIM, params);
    catch ME
        IM([1 end],:,:) = nan;
        IM(:,[1 end],:) = nan;
        [IMc, peaks, actAcc] = localizeSources_vIM(IM, vIM, params);
    end
end