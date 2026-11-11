function s = loadStructFromH5(filename, skipDatasets)
%LOADSTRUCTFROMH5 Read an HDF5 file written by SAVESTRUCTTOH5 into a struct.
%   s = LOADSTRUCTFROMH5(filename) mirrors the hierarchy written by
%   SAVESTRUCTTOH5. HDF5 groups become nested scalar structs; string
%   datasets become cell arrays of char, except scalar path strings
%   (e.g. datadr, savedr) which become char. filename and output path
%   grids stay cell with the same shape as the dataset (including 1x1).
%
%   s = LOADSTRUCTFROMH5(filename, skipDatasets) does not read the datasets
%   whose full paths are listed in skipDatasets (char or cellstr), e.g.
%   '/slap2/varFacDS', the multi-GB variance movie of an alignment file, when
%   only its metadata are needed. Skipped datasets are absent from s.

if nargin < 1 || isempty(filename)
    error('loadStructFromH5:MissingFilename', 'A source filename must be provided.');
end
if ~exist(filename, 'file')
    error('loadStructFromH5:FileNotFound', 'File not found: %s', filename);
end
if nargin < 2
    skipDatasets = {};
end
skipDatasets = cellstr(skipDatasets);

info = h5info(filename);
s = readGroup(filename, info, skipDatasets);
if isfield(s, 'slap2_info') && isstruct(s.slap2_info) && isfield(s.slap2_info, 'ref_stack')
    rs = s.slap2_info.ref_stack;
    for fn = fieldnames(rs)'
        tok = regexp(fn{1}, '^DMD(\d+)$', 'tokens', 'once');
        if ~isempty(tok)
            pathFn = ['Path' tok{1}];
            if ~isfield(rs, pathFn)
                rs.(pathFn) = rs.(fn{1});
            end
            rs = rmfield(rs, fn{1});
        end
    end
    s.slap2_info.ref_stack = rs;
end
if isfield(s, 'slap2') && isstruct(s.slap2)
    if isfield(s.slap2, 'Z') && ~isfield(s.slap2, 'Z_depths')
        s.slap2.Z_depths = s.slap2.Z;
        s.slap2 = rmfield(s.slap2, 'Z');
    end
end
end


function s = readGroup(filename, grpInfo, skipDatasets)
%READGROUP Recursively convert an h5info group into a scalar struct.
s = struct();

for ix = 1:numel(grpInfo.Datasets)
    ds = grpInfo.Datasets(ix);
    if strcmp(grpInfo.Name, '/')
        path = ['/' ds.Name];
    else
        path = [grpInfo.Name '/' ds.Name];
    end
    if any(strcmp(path, skipDatasets))
        continue
    end
    val = h5read(filename, path);
    if isstring(val)
        if h5TrialTableFilenameCellField(ds.Name)
            val = reshape(cellstr(val(:)), size(val));
        elseif isscalar(val)
            val = char(val);
        else
            val = cellstr(val);
        end
    elseif ischar(val) && h5TrialTableFilenameCellField(ds.Name) && isvector(val)
        % Rare: string dataset read as char (e.g. single-trial path)
        val = {strtrim(val(:)')};
    end
    s.(ds.Name) = val;
end

for ix = 1:numel(grpInfo.Groups)
    grp = grpInfo.Groups(ix);
    parts = strsplit(grp.Name, '/');
    grpName = parts{end};
    s.(grpName) = readGroup(filename, grp, skipDatasets);
end
end

function tf = h5TrialTableFilenameCellField(name)
persistent known
if isempty(known)
    known = {'filename', 'fn_reg_ds', 'fn_raw', 'fn_adata', 'fn_ann'};
end
tf = any(strcmp(name, known));
end
