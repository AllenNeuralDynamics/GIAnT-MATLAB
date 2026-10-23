function params = StripRegistration(fullPathToTrialTable, paramsIn)
%STRIPREGISTRATION Strip alignment using a pre-built trial_table.h5.
%   params = STRIPREGISTRATION() prompts for trial_table.h5 (see BUILDTRIALTABLE).
%   params = STRIPREGISTRATION(fullPathToTrialTable)
%   params = STRIPREGISTRATION(fullPathToTrialTable, paramsIn)  % same pattern as MultiRoiRegistration

if ~nargin || isempty(fullPathToTrialTable)
    [fn, trialtabledr] = uigetfile('*.h5', 'Select a trial_table file', '*trial_table*.h5');
    fullPathToTrialTable = fullfile(trialtabledr,fn);
    if isequal(fn, 0)
        error('StripRegistration:NoTrialTable', 'A trial_table.h5 file is required. Build one with buildTrialTable first.');
    end
else
    [trialtabledr, fn, ext] = fileparts(fullPathToTrialTable);
    fn = [fn ext];
end

trialTable = loadStructFromH5(fullfile(trialtabledr, fn));

mocodr = fullfile(trialTable.savedr, 'motion_correction');
if ~exist(mocodr, 'dir')
    mkdir(mocodr);
end
dataDr = trialTable.datadr;

%PARAMETER SETTING
if nargin > 1
    if ischar(paramsIn) || isstring(paramsIn)
        paramsIn = jsondecode(char(paramsIn));
    end
    params = setParams('StripRegistration', paramsIn);
else
    params = setParams('StripRegistration');
end

%set up parallelization
if params.nWorkers>1
    p = gcp('nocreate'); % If no pool, do not create new one.
    if isempty(p)
        poolsize = 0;
    else
        poolsize = p.NumWorkers;
    end
    nWorkers = min(params.nWorkers, numel(trialTable.filename));
    if poolsize<nWorkers
        delete(gcp('nocreate'));
        if params.nWorkers<4
            warning('You are using few parallel workers!');
        end
        parpool('processes',nWorkers); %limit the number of workers to avoid running out of RAM
    end
end

%align in parallel
nTrials= numel(trialTable.filename);
fnRegDS = cell(1, nTrials); fnRaw = cell(1, nTrials); fnAdata = cell(1, nTrials);
regFail = false(1, nTrials);
if params.nWorkers>1
    parfor f_ix = 1:nTrials
        [fnRegDS{f_ix}, fnRaw{f_ix}, fnAdata{f_ix}, regFail(f_ix)] = alignAsync(dataDr, mocodr, trialTable, params, f_ix);
    end
else
    for f_ix = 1:nTrials
        [fnRegDS{f_ix}, fnRaw{f_ix}, fnAdata{f_ix}, regFail(f_ix)] = alignAsync(dataDr, mocodr, trialTable, params, f_ix);
    end
end

if ~isfield(trialTable, 'motion_correction') || isempty(trialTable.motion_correction)
    trialTable.motion_correction = struct();
end
trialTable.motion_correction.fn_reg_ds = fnRegDS;
trialTable.motion_correction.fn_raw = fnRaw;
trialTable.motion_correction.fn_adata = fnAdata;
trialTable.motion_correction.registration_failed = regFail;
trialTable.motion_correction.align_params = params;
saveStructToH5(trialTable, fullPathToTrialTable);

disp('done StripRegistration.')

end


function [fnDS, fnRaw, fnAdata, registrationFailed] = alignAsync(dataDr, outDr, trialTable, params, f_ix)
maxshift = params.maxshift;
dsFac = 2^params.ds_time; params.dsFac = dsFac;
registrationFailed = false;

fn = trialTable.filename{f_ix};
[~,fn, ext] = fileparts(fn);
fnstem = fn;
fn = strcat(fn,ext);

disp(['Aligning: ' [dataDr filesep fn]])

if endsWith(fn, '.h5')
    desc = h5info([dataDr filesep fn]);
    Ad = single(readStack([dataDr filesep fn], [], 'Dataset', desc.Datasets(1).Name));
else
    [Ad, meta, desc] = ScanImageTiffWrapper([dataDr filesep fn]);
    Ad = single(Ad);
end
try
    evalc(desc{1});

    metaLines = strsplit(meta, '\n');
    for lineIx = 1:length(metaLines)
        try
            evalc([metaLines{lineIx} ';']);
        catch
        end
    end
    numChannels = length(SI.hChannels.channelSave);
    selCh = 1:numChannels;
catch  % assume 1 channel if there is no metadata (e.g. simulated data)
    numChannels = 1;
    selCh = 1;
end

try
    pat = "frameTimestamps_sec = " + digitsPattern + "." + digitsPattern;
    for frame = 1:10*numChannels %compute the framerate from the metadata by reading a few frames
        E = extract(desc{frame}, pat);
        timestamp(frame) = str2double(E{1}(23:end)); %#ok<AGROW>
    end
    frametime = median(diff(timestamp(1:numChannels:end)));
catch  % use default frametime if there is no metadata (e.g. simulated data)
    if params.frameRate > 0
        frametime = 1/params.frameRate;
        fprintf('Warning! Failed to compute framerate from metadata. Using user input frametime=%0.4f\n',frametime);
    else
        frametime = 0.0023;
        disp('Warning! Failed to compute framerate from metadata. Using default frametime=0.0023')
    end
end
aData.numChannels = numChannels;
aData.frametime = frametime;
aData.alignHz = 1/frametime/dsFac;

Ad = permute(reshape(Ad, size(Ad,1), size(Ad,2), numChannels, []), [2 1 3 4]);
Ad = Ad(params.removeLines+1:end,:,:,:);

%subtract baseline
tsamp = round(linspace(1,size(Ad,4), 500));
BG = min(1000,prctile(Ad(:,:,:,tsamp), 10, [1 2 4]));
Ad = Ad-min(BG,1000);

%make an initial template with normcorr
initFrames = 1000;
framesToRead = initFrames * dsFac;
try
    Y = downsampleTime(Ad(:,:,:,1:framesToRead), params.ds_time);
catch
    disp(['Your file was too short:' fn])
    fnDS = []; fnRaw=[]; fnAdata = [];
    registrationFailed = true;
    return
end
sz = size(Ad);
Yhp = squeeze(sum(reshape(Y(:,:,selCh,:),size(Y,1),size(Y,2),[],size(Y,4)),3));

% find the most correlated ~100 frames within the first initFrames to
% use as an initial template
rho = corr(reshape(Yhp,[],size(Yhp,3)));
dist_matrix = 1 - rho;
Z = linkage(squareform(dist_matrix),'average');

cutoff = 0.01;
min_cluster_size = 100;
clusters = [];
max_cutoff = 2.0;
while isempty(clusters) || all(cellfun(@numel, clusters) < min_cluster_size)
    cutoff = cutoff + 0.01;
    if cutoff > max_cutoff
        error('Could not find a cluster with at least %d samples',min_cluster_size);
    end
    T = cluster(Z, 'cutoff', cutoff, 'criterion', 'distance');
    clusters = arrayfun(@(x) find(T == x), unique(T), 'UniformOutput', false);
end

max_mean_corr = -Inf;
for i = 1:length(clusters)
    cluster_indices = clusters{i};
    if numel(cluster_indices) >= min_cluster_size
        mean_corr = mean(mean(rho(cluster_indices,cluster_indices)));
        if mean_corr > max_mean_corr
            max_mean_corr = mean_corr;
            best_cluster = cluster_indices;
        end
    end
end

options_rigid = NoRMCorreSetParms('d1',size(Y,1),'d2',size(Y,2),'bin_width',initFrames,'max_shift',maxshift,'us_fac',4,'init_batch',initFrames, 'correct_bidir', false);
F = normcorre(Yhp,options_rigid,mean(Yhp(:,:,best_cluster),3));
F = mean(F,3); %fixed image
template = nan(2*maxshift+sz(1), 2*maxshift+sz(2));
templateFull = template;
T0 = template; T00 = zeros(size(template)); templateCt = zeros(size(template));
T0(maxshift+(1:sz(1)), maxshift+(1:sz(2)))=F;
clear Y Yhp;

initR = 0; initC = 0;
nDSframes= floor(sz(4)./(dsFac)); %number of downsampled frames
motionDSr = nan(1,nDSframes); motionDSc = nan(1,nDSframes); %matrices to store the inferred motion
[viewR, viewC] = ndgrid((1:(sz(1)+2*maxshift))-maxshift, (1:(sz(2)+2*maxshift))-maxshift); %view matrices for interpolation

%for spatial downsampling, used to calculate alignment quality
dsTimes = 2;
dsSz = floor(size(template)./(2^dsTimes));
A_ds = nan([dsSz nDSframes],'single');

disp('Registering:');
for DSframe = 1:nDSframes
    readFrames = (DSframe-1)*(dsFac) + (1:(dsFac));

    M = downsampleTime(Ad(:,:,:, readFrames), params.ds_time);
    M = squeeze(sum(reshape(M(:,:,selCh,:),size(M,1),size(M,2),[],size(M,4)),3)); %merge colors

    if ~mod(DSframe, 1000)
        disp([int2str(DSframe) ' of ' int2str(nDSframes)]);
    end

    Ttmp = mean(cat(3, T0, T00,template),3, 'omitnan');
    T = Ttmp(maxshift-initR + (1:sz(1)), maxshift-initC+(1:sz(2)));

    output = dftregistration_clipped(fft2(M),fft2(single(T)),4, params.clipShift);
    motionDSr(DSframe) = initR+output(3);
    motionDSc(DSframe) = initC+output(4);

    if sqrt((motionDSr(DSframe)/sz(1)).^2 + (motionDSc(DSframe)/sz(2)).^2) > 0.75^2
        Mfull = interp2(1:sz(2), 1:sz(1), M,viewC, viewR, 'linear', nan);
        [motion, ~] = xcorr2_nans(Mfull,Ttmp,[initR;initC],50);

        motionDSr(DSframe) = motion(1);
        motionDSc(DSframe) = motion(2);
    end

    if abs(motionDSr(DSframe))<maxshift && abs(motionDSc(DSframe))<maxshift
        A = interp2(1:sz(2), 1:sz(1), M,viewC+motionDSc(DSframe), viewR+motionDSr(DSframe), 'linear', nan);

        %downsample in space for alignment quality stats
        dsTmp = A;
        for dsIx = 1:dsTimes
            dsTmp = dsTmp(1:2:2*floor(end/2), 1:2:2*floor(end/2)) + dsTmp(1:2:2*floor(end/2), 2:2:2*floor(end/2)) + dsTmp(2:2:2*floor(end/2), 1:2:2*floor(end/2)) + dsTmp(2:2:2*floor(end/2), 2:2:2*floor(end/2));
        end
        A_ds(:,:,DSframe) = dsTmp;

        templateFull = sum(cat(3,templateFull .* templateCt, A),3,"omitnan");
        templateCt = templateCt + ~isnan(A);
        templateFull = templateFull ./ templateCt;
        template = templateFull;
        template(templateCt < 100) = nan;

        initR = round(motionDSr(DSframe));
        initC = round(motionDSc(DSframe));
    else
        motionDSr(DSframe) = initR;
        motionDSc(DSframe) = initC;
    end
end

%compute alignment error stats
offset = max(0, prctile(template(:), 25));
recNegErr = nan(1,nDSframes);
nChunks = ceil(nDSframes./(aData.alignHz*10)); %align 10 seconds at a time
chunkEdges = round(linspace(1, nDSframes+1, nChunks));
for chunkIx = 1:length(chunkEdges)-1
    t_ixs = chunkEdges(chunkIx):chunkEdges(chunkIx+1)-1;
    chunkData = A_ds(:,:,t_ixs);
    template = median(chunkData, 3, 'omitmissing');
    nanFrac = mean(isnan(chunkData), 3);
    template(nanFrac>0.2) = nan;
    template_gamma = sqrt(max(0,template)+offset);
    % chunkData*0: broadcast 2D template to 3D without repmat; nan where frame is nan
    recNegErr(1,t_ixs) = sqrt(squeeze(mean(max(0, (template-chunkData)./template_gamma).^2, [1 2], 'omitnan')./mean(max(0, (template./template_gamma).^2 + chunkData * 0, 'includenan'), [1 2], 'omitnan')));
end

%upsample the shifts and compute a tighter field of view
tDS = (1:nDSframes).*(dsFac) - (2^(params.ds_time-1)) + 0.5;
motionC = interp1(tDS, motionDSc, 1:((2^params.ds_time)*nDSframes), 'pchip','extrap'); %upsample the movement to the original timebase
motionR = interp1(tDS, motionDSr, 1:((2^params.ds_time)*nDSframes), 'pchip', 'extrap'); %upsample the movement to the original timebase
maxshiftC = max(abs(motionC)); maxshiftR = max(abs(motionR));
[viewR, viewC] = ndgrid((1:(sz(1)+2*maxshiftR))-maxshiftR, (1:(sz(2)+2*maxshiftC))-maxshiftC); %view matrices for interpolation

pixelscale = 4e4; %PIXEL SIZE IN DOTS PER CM

%save a downsampled aligned recording
fnDS = [fnstem '_REGISTERED_DOWNSAMPLED-' int2str(dsFac) 'x.tif'];
if ~params.saveTif
    fnDS = strrep(fnDS, '.tif', '.h5');
end
Bsum = zeros([size(viewR') numChannels]);
Bcount = zeros([size(viewR') numChannels]);

tiffSave = zeros([size(viewR, [2 1]) nDSframes*numChannels], 'single');

for DSframe = 1:nDSframes
    readFrames = (DSframe-1)*(dsFac) + (1:(dsFac));
    YY = downsampleTime(Ad(:,:,:, readFrames), params.ds_time);
    for ch = 1:numChannels
        B =  interp2(1:sz(2), 1:sz(1), YY(:,:,ch),viewC+motionDSc(DSframe), viewR+motionDSr(DSframe), 'linear', nan)';
        tiffSave(:,:,(DSframe-1)*numChannels+ch) = single(B);

        Bcount(:,:,ch) = Bcount(:,:,ch) + ~isnan(B);
        B(isnan(B)) = 0;
        Bsum(:,:,ch) = Bsum(:,:,ch)+double(B);
    end
end

nanRows = mean(mean(isnan(tiffSave),3),2) == 1;
nanCols = mean(mean(isnan(tiffSave),3),1) == 1;

tiffSave(nanRows,:,:) = [];
tiffSave(:,nanCols,:) = [];

[h, w, nt] = size(tiffSave);
if nt ~= nDSframes * numChannels
    error('StripRegistration:stackSize', 'Expected %d pages (%d ch x %d frames), got %d.', nDSframes * numChannels, numChannels, nDSframes, nt);
end
% Per-channel mean over downsampled aligned frames, numChannels x H x W
meanIM = permute(single(mean(reshape(tiffSave, h, w, numChannels, nDSframes), 4, 'omitnan')), [3 1 2]);

if params.saveTif
    writeStackBigTiff(tiffSave, [outDr filesep fnDS], pixelscale);
else
    h5fn = [outDr filesep strrep(fnDS,'.tif','.h5')];
    h5sz = size(tiffSave);
    h5create(h5fn, '/data', h5sz, 'Datatype', 'single', 'Deflate', 4, 'ChunkSize', [h5sz(1) h5sz(2) min(500, h5sz(3))]);
    h5write(h5fn, '/data', tiffSave);
end
clear('tiffSave');

%save an average image for each channel
for ch = 1:numChannels
    Bmean = Bsum(:,:,ch)./Bcount(:,:,ch);
    minV = prctile(Bmean(~isnan(Bmean(:))), 10);
    maxV = prctile(Bmean(~isnan(Bmean(:))), 99.9);
    Bmean = uint8(255*sqrt(max(0,(Bmean-minV)./(maxV-minV))));
    Bmean(nanRows,:) = [];
    Bmean(:,nanCols) = [];
    fnwrite = [fnstem '_REGISTERED_AVG_CH' num2str(ch) '_8bit.tif'];
    writeStackBigTiff(single(Bmean), [outDr filesep fnwrite], pixelscale);
end

fnRaw = [fnstem '_REGISTERED_RAW.tif'];
if ~params.saveTif
    fnRaw = strrep(fnRaw, '.tif', '.h5');
end
tiffSave = zeros([size(viewR, [2 1]) - [sum(nanRows) sum(nanCols)] length(motionC)*numChannels], 'single');
for frame = 1:length(motionC)
    for ch = 1:numChannels
        B = interp2(1:sz(2), 1:sz(1), Ad(:,:,ch,frame),viewC+motionC(frame), viewR+motionR(frame), 'linear', nan)';
        B(nanRows,:) = [];
        B(:,nanCols) = [];
        tiffSave(:,:,(frame-1)*numChannels+ch) = single(B);
    end
end

if params.saveTif
    writeStackBigTiff(single(tiffSave), [outDr filesep fnRaw], pixelscale);
else
    h5fn = [outDr filesep strrep(fnRaw,'.tif','.h5')];
    h5sz = size(tiffSave);
    h5create(h5fn, '/data', h5sz, 'Datatype', 'single', 'Deflate', 4, 'ChunkSize', [h5sz(1) h5sz(2) min(500, h5sz(3))]);
    h5write(h5fn, '/data', single(tiffSave));
end
clear('tiffSave')

%save alignment metadata
fnAdata = [fnstem '_ALIGNMENTDATA.h5'];
toSave = struct();
toSave.numChannels = aData.numChannels;
toSave.meanIM = meanIM;
toSave.frametime = aData.frametime;
toSave.alignHz = aData.alignHz;
toSave.motionC = motionC;
toSave.motionR = motionR;
toSave.motionDSc = motionDSc;
toSave.motionDSr = motionDSr;
toSave.recNegErr = recNegErr;
saveStructToH5(toSave, [outDr filesep fnAdata]);
end


function Y = downsampleTime(Y, ds_time)
for ix = 1:ds_time
    Y = Y(:,:,:,1:2:(2*floor(end/2)))+ Y(:,:,:,2:2:end);
end
Y = Y./2.^ds_time;
end

function writeStackBigTiff(img, filepath, pixelscale)
fTIF = Fast_BigTiff_Write(filepath, pixelscale, 0);
if ndims(img) <= 2
    fTIF.WriteIMG(img);
else
    for k = 1:size(img, 3)
        fTIF.WriteIMG(img(:, :, k));
    end
end
fTIF.close;
end