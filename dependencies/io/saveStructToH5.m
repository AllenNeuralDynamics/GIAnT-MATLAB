function saveStructToH5(s, filename)
%SAVESTRUCTTOH5 Write a (possibly nested) scalar struct to an HDF5 file.
%   SAVESTRUCTTOH5(s, filename) recursively maps MATLAB struct fields onto
%   the HDF5 hierarchy:
%     - Scalar struct fields become HDF5 groups.
%     - Cell arrays of char rows/string scalars become variable-length
%       string datasets (empty elements as "").
%     - char rows and string scalars/arrays become string datasets.
%     - Numeric and logical arrays become numeric datasets (logicals as int8).
%     - Empty fields are skipped; a group exists only above a dataset.
%
%   The destination file is overwritten if it already exists; HDF5 datasets
%   cannot be redefined in place so a fresh file is required.
%
%   Dataset sizes match MATLAB size(val). The file includes /row_major
%   (uint8): 0 = column-major (MATLAB layout). See README.md
%   "Reading H5 outside MATLAB".
%
%   The whole file is written through one HDF5 file handle. h5create and
%   h5write reopen the file for every dataset, which on a network share
%   made a trial table (72 datasets) take 5-6 s instead of under 1 s. Every
%   value is checked before the file is touched, so an unsupported value
%   leaves an existing file as it was. Creating the file and one read-back
%   check of the finished file are retried briefly, for transient errors of
%   network storage. Files are the same as the h5create/h5write writer's:
%   paths, sizes, native types, variable-length UTF-8 strings, contiguous
%   layout.

if nargin < 2 || isempty(filename)
    error('saveStructToH5:MissingFilename', 'A destination filename must be provided.');
end
if ~isstruct(s) || ~isscalar(s)
    error('saveStructToH5:UnsupportedInput', 'Input must be a scalar struct.');
end

% /row_major is written below; a struct read with loadStructFromH5 (e.g. a
% trial table being updated) also carries it as a field
if isfield(s, 'row_major')
    s = rmfield(s, 'row_major');
end
validateGroup('', s);

fid = createFile(filename);
lcpl = H5P.create('H5P_LINK_CREATE');
closer = onCleanup(@() closeIds(fid, lcpl));
H5P.set_create_intermediate_group(lcpl, 1); %groups are made with their first dataset, as h5create does

written = writeGroup(fid, lcpl, '', s, cell(0, 2)); %{path, size} per dataset
written = writeDataset(fid, lcpl, '/row_major', uint8(0), written);
H5P.close(lcpl);
H5F.close(fid); %not in the cleanup, so that a failed flush is reported
clear closer

verifyReadBack(filename, written);
end


%% validation (no file I/O)

function validateGroup(basePath, s)
fields = fieldnames(s);
for ix = 1:numel(fields)
    validateValue([basePath '/' fields{ix}], s.(fields{ix}));
end
end


function validateValue(path, val)
if isstruct(val)
    if ~isscalar(val)
        error('saveStructToH5:NonScalarStruct', ...
            'Non-scalar structs are not supported (at %s).', path);
    end
    validateGroup(path, val);
elseif iscell(val)
    % each element becomes one string of the dataset
    if ~all(cellfun(@isTextScalar, val(:)))
        error('saveStructToH5:UnsupportedCell', ...
            'Only cell arrays of char rows or string scalars are supported (at %s).', path);
    end
elseif ischar(val)
    if ~isempty(val) && ~isrow(val)
        error('saveStructToH5:UnsupportedChar', ...
            'Only char row vectors are supported (at %s); use a string array.', path);
    end
elseif isnumeric(val)
    if isempty(nativeTypeName(class(val)))
        error('saveStructToH5:UnsupportedNumeric', ...
            'Unsupported numeric class ''%s'' at %s.', class(val), path);
    end
    if issparse(val) || ~isreal(val)
        error('saveStructToH5:UnsupportedNumeric', ...
            'Sparse and complex arrays are not supported (at %s).', path);
    end
elseif islogical(val)
    if issparse(val) %written as int8, which cannot be sparse
        error('saveStructToH5:UnsupportedLogical', ...
            'Sparse logical arrays are not supported (at %s).', path);
    end
elseif ~isstring(val)
    error('saveStructToH5:UnsupportedType', ...
        'Unsupported value type ''%s'' at %s.', class(val), path);
end
end


function tf = isTextScalar(x)
%ISTEXTSCALAR A cell element that writes as one string ("" if empty)
tf = isempty(x) || (ischar(x) && isrow(x)) || (isstring(x) && isscalar(x));
end


%% writing

function written = writeGroup(fid, lcpl, basePath, s, written)
fields = fieldnames(s);
for ix = 1:numel(fields)
    written = writeValue(fid, lcpl, [basePath '/' fields{ix}], s.(fields{ix}), written);
end
end


function written = writeValue(fid, lcpl, path, val, written)
if isstruct(val)
    written = writeGroup(fid, lcpl, path, val, written);
elseif iscell(val)
    if isempty(val)
        return
    end
    sval = strings(size(val));
    for k = 1:numel(val)
        if ~isempty(val{k})
            sval(k) = string(val{k});
        end
    end
    written = writeDataset(fid, lcpl, path, sval, written);
elseif ischar(val)
    written = writeDataset(fid, lcpl, path, string(val), written); %[1 1]
elseif isstring(val)
    written = writeDataset(fid, lcpl, path, val, written);
elseif islogical(val)
    if ~isempty(val)
        written = writeDataset(fid, lcpl, path, int8(val), written);
    end
elseif ~isempty(val) %numeric
    val = gather(val); %collect gpuArray / codistributed
    % Round-trip through cast to materialize a plain array of the declared
    % class without applying per-element arithmetic.
    written = writeDataset(fid, lcpl, path, cast(val, class(val)), written);
end
end


function written = writeDataset(fid, lcpl, path, val, written)
%WRITEDATASET Create and write one dataset of size(val) on the open file.
% HDF5 dims are MATLAB's reversed, as with h5create/h5write.
sz = size(val);
if isstring(val)
    tid = H5T.copy('H5T_C_S1');
    H5T.set_size(tid, 'H5T_VARIABLE');
    H5T.set_cset(tid, 'H5T_CSET_UTF8');
    data = cellstr(val);
else
    tid = H5T.copy(nativeTypeName(class(val)));
    data = val;
end
sid = H5S.create_simple(numel(sz), fliplr(sz), []);
did = H5D.create(fid, path, tid, sid, lcpl, 'H5P_DEFAULT', 'H5P_DEFAULT');
H5D.write(did, tid, 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', data);
H5D.close(did);
H5S.close(sid);
H5T.close(tid);
written(end+1, :) = {path, sz};
end


function name = nativeTypeName(cls)
switch cls
    case 'double', name = 'H5T_NATIVE_DOUBLE';
    case 'single', name = 'H5T_NATIVE_FLOAT';
    case 'int8',   name = 'H5T_NATIVE_INT8';
    case 'uint8',  name = 'H5T_NATIVE_UINT8';
    case 'int16',  name = 'H5T_NATIVE_INT16';
    case 'uint16', name = 'H5T_NATIVE_UINT16';
    case 'int32',  name = 'H5T_NATIVE_INT32';
    case 'uint32', name = 'H5T_NATIVE_UINT32';
    case 'int64',  name = 'H5T_NATIVE_INT64';
    case 'uint64', name = 'H5T_NATIVE_UINT64';
    otherwise,     name = '';
end
end


%% file create / close / verify (retried)

function fid = createFile(filename)
%CREATEFILE Create (truncate) the file, retrying transient failures.
% Truncating replaces delete-then-create, which can race a handle left on
% the deleted file (virus scanner, sync client). H5F_CLOSE_STRONG closes
% any objects still open when the file is closed.
maxAttempts = 6;
for attempt = 1:maxAttempts
    fapl = H5P.create('H5P_FILE_ACCESS');
    try
        H5P.set_fclose_degree(fapl, 'H5F_CLOSE_STRONG');
        fid = H5F.create(filename, 'H5F_ACC_TRUNC', 'H5P_DEFAULT', fapl);
        H5P.close(fapl);
        return
    catch ME
        H5P.close(fapl);
        if attempt == maxAttempts
            error('saveStructToH5:CreateFailed', ...
                'Could not create %s after %d attempts: %s', filename, maxAttempts, ME.message);
        end
        pause(backoffSeconds(attempt));
    end
end
end


function closeIds(fid, lcpl)
%CLOSEIDS onCleanup after an error: close the file (and with it any open
%objects) and the link property list
try
    if H5I.is_valid(fid)
        H5F.close(fid);
    end
    if H5I.is_valid(lcpl)
        H5P.close(lcpl);
    end
catch
end
end


function verifyReadBack(filename, written)
%VERIFYREADBACK Reopen the finished file once and check every dataset.
% A failure to reopen is retried; a missing or misshapen dataset is a bug,
% not a storage delay, and is reported at once.
maxChecks = 10;
for checkIx = 1:maxChecks
    try
        found = collectDatasets(h5info(filename), containers.Map('KeyType', 'char', 'ValueType', 'any'));
    catch ME
        if checkIx == maxChecks
            error('saveStructToH5:VerifyFailed', ...
                'Wrote %s but could not read it back after %d attempts: %s', ...
                filename, maxChecks, ME.message);
        end
        pause(backoffSeconds(checkIx));
        continue
    end
    for ix = 1:size(written, 1)
        path = written{ix, 1};
        if ~isKey(found, path)
            error('saveStructToH5:DatasetMissing', ...
                'Dataset %s is missing from %s after writing.', path, filename);
        end
        if ~isequal(double(found(path)), double(written{ix, 2}))
            error('saveStructToH5:DatasetSizeMismatch', ...
                'Dataset %s has size [%s], expected [%s].', path, ...
                num2str(found(path)), num2str(written{ix, 2}));
        end
    end
    return
end
end


function found = collectDatasets(grpInfo, found)
%COLLECTDATASETS Map of dataset path -> Dataspace.Size for an h5info tree.
for ix = 1:numel(grpInfo.Datasets)
    if strcmp(grpInfo.Name, '/')
        path = ['/' grpInfo.Datasets(ix).Name];
    else
        path = [grpInfo.Name '/' grpInfo.Datasets(ix).Name];
    end
    found(path) = grpInfo.Datasets(ix).Dataspace.Size;
end
for ix = 1:numel(grpInfo.Groups)
    found = collectDatasets(grpInfo.Groups(ix), found);
end
end


function s = backoffSeconds(attempt)
s = min(0.8, 0.05 * 2^(attempt - 1)); %0.05, 0.1, 0.2, 0.4, 0.8 s
end
