function [IM, fn] = getAnnotationImage(trialTable, keepTrials, mocodr, DMDix)
%GETANNOTATIONIMAGE Mean image of the longest valid trial, the space ROIs are annotated in.
%
% Returns IM (H x W) and the registered file stem fn, or empty values if the
% imaging path has no valid trials. The longest trial (the most registered
% movie data) gives the cleanest mean image; short trials, often the first
% one, can be too noisy to annotate on or to register post-acquisition
% annotations to. Ties go to the first trial.
% Prefers meanIM from the alignment H5 and falls back to averaging the
% registered movie.

IM = [];
fn = '';
validTrials = find(keepTrials(DMDix,:));
if isempty(validTrials)
    return
end
nBytes = zeros(size(validTrials));
for ix = 1:numel(validTrials)
    nBytes(ix) = movieDataBytes(fullfile(mocodr, ...
        trialTable.motion_correction.fn_reg_ds{DMDix, validTrials(ix)}));
end
[~, best] = max(nBytes);
annTrial = validTrials(best);

[~, fn, ext] = fileparts(trialTable.motion_correction.fn_reg_ds{DMDix,annTrial});
adataFn = trialTable.motion_correction.fn_adata{DMDix, annTrial};
aDataAnnot = loadStructFromH5(fullfile(mocodr, adataFn), '/slap2/varFacDS');
if isfield(aDataAnnot, 'meanIM') && ~isempty(aDataAnnot.meanIM)
    IM = squeeze(mean(aDataAnnot.meanIM, 1, 'omitnan')); % numChannels x H x W -> H x W
elseif strcmp(ext, '.h5')
    movieFn = fullfile(mocodr, [fn ext]);
    desc = h5info(movieFn);
    IMh5 = h5read(movieFn, ['/', desc.Datasets(1).Name]); % H x W x (channels x frames)
    IM = mean(IMh5, 3, 'omitnan');
else
    [IMtif, ~] = ScanImageTiffWrapper(fullfile(mocodr, [fn ext]));
    IM = squeeze(mean(IMtif, [3 4], 'omitnan'));
end
end

function nBytes = movieDataBytes(movieFn)
% Size of a registered movie's data, which grows with its frames. h5 movies
% may be compressed, and their file size then depends on content as well,
% so they are sized by their dataset (single); other movies by their file.
% Missing files count as 0.
nBytes = 0;
d = dir(movieFn);
if isempty(d)
    return
end
nBytes = d(1).bytes;
if endsWith(movieFn, '.h5')
    try
        desc = h5info(movieFn);
        nBytes = 4*prod(desc.Datasets(1).Dataspace.Size);
    catch
        % unreadable: keep the file size
    end
end
end
