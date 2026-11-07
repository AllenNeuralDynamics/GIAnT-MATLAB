function [trialTable, keepTrials] = verifyFiles(fn,dr)

trialTable = loadStructFromH5([dr filesep fn]);

mocodr = [trialTable.savedr filesep 'motion_correction'];
datadr = trialTable.datadr;

nTrials = size(trialTable.true_trial_ix,2);
nDMDs = size(trialTable.filename,1);
keepTrials = true(nDMDs, nTrials);
useRegFailTable = isfield(trialTable, 'motion_correction') ...
    && isfield(trialTable.motion_correction, 'registration_failed') ...
    && size(trialTable.motion_correction.registration_failed, 1) == nDMDs ...
    && size(trialTable.motion_correction.registration_failed, 2) == nTrials;

%populate source_extraction.fn_raw with the file that source extraction
%should read raw data from (.dat for SLAP2, registered-raw tif for non-SLAP2)
trialTable.source_extraction.fn_raw = cell(nDMDs, nTrials);

for trialIx = nTrials:-1:1
    for DMDix = 1:nDMDs
        [~, tiffFn, ext] = fileparts(trialTable.motion_correction.fn_reg_ds{DMDix,trialIx});
        if isempty(dir([mocodr filesep tiffFn '.tif'])) & isempty(dir([mocodr filesep tiffFn '.h5']))
            disp(['Missing tiff or h5 file:' tiffFn])
            keepTrials(DMDix,trialIx) = false;
        else
            trialTable.motion_correction.fn_reg_ds{DMDix,trialIx} = [tiffFn, ext];
        end

        [~, alignFn] = fileparts(trialTable.motion_correction.fn_adata{DMDix,trialIx});
        if ~exist([mocodr filesep alignFn '.h5'], 'file')
            disp(['Missing alignData file:' alignFn])
            keepTrials(DMDix,trialIx) = false;
        else
            if useRegFailTable
                regFailed = trialTable.motion_correction.registration_failed(DMDix, trialIx);
            else
                aData = loadStructFromH5([mocodr filesep alignFn '.h5'], '/slap2/varFacDS');
                regFailed = isfield(aData, 'registrationFailed') && aData.registrationFailed;
                clear aData
            end
            if regFailed
                disp(['Registration failed for file:' alignFn])
                keepTrials(DMDix,trialIx) = false;
            else
                trialTable.motion_correction.fn_adata{DMDix,trialIx} = [alignFn '.h5'];
            end
        end

        sourceFn = trialTable.filename{DMDix,trialIx};
        if ~exist([datadr filesep sourceFn], 'file')
            disp(['Missing source data file:' sourceFn])
            keepTrials(DMDix, trialIx) = false;
        end
        [~,~,sourceExt] = fileparts(sourceFn);
        if strcmpi(sourceExt, '.dat')
            %SLAP2: raw data for extraction is the source .dat file in datadr
            trialTable.source_extraction.fn_raw{DMDix,trialIx} = sourceFn;
        else
            %non-SLAP2: raw data for extraction is the registered-raw file
            %produced by StripRegistration and stored in motion_correction
            rawFn = trialTable.motion_correction.fn_raw{DMDix,trialIx};
            if ~exist([mocodr filesep rawFn], 'file')
                disp(['Missing raw data file:' rawFn])
                keepTrials(DMDix, trialIx) = false;
            end
            trialTable.source_extraction.fn_raw{DMDix,trialIx} = rawFn;
        end
    end
end
if ~all(keepTrials)
    disp(['Files were missing for ' int2str(sum(~keepTrials(:))) ' recordings; likely failed alignments. Proceeding without them.']);
end
if ~any(keepTrials)
    error('All trials were rejected due to missing alignment files!');
end
end
