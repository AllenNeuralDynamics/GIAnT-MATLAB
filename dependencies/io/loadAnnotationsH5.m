function [ROIs, schemaVersion] = loadAnnotationsH5(filename)
%LOADANNOTATIONSH5 Load manual ROI annotations from an HDF5 file.
%
% schemaVersion is the file's /schema_version (1 for unversioned files);
% see annotationsSchemaVersion.

if nargin < 1 || isempty(filename)
    error('loadAnnotationsH5:MissingFilename', 'A source filename must be provided.');
end
if ~exist(filename, 'file')
    error('loadAnnotationsH5:FileNotFound', 'File not found: %s', filename);
end

schemaVersion = double(readNumericOrDefault(filename, '/schema_version', uint32(1)));
if schemaVersion > annotationsSchemaVersion()
    warning('loadAnnotationsH5:NewerSchema', ...
        ['%s uses annotation schema version %d, newer than this GIAnT supports (%d); ' ...
         'update GIAnT-MATLAB. Fields it does not know about are ignored.'], ...
        filename, schemaVersion, annotationsSchemaVersion());
end
coordsZeroIndexed = logical(readNumericOrDefault(filename, '/coords_zero_indexed', uint8(0)));

info = h5info(filename);
pathNames = {};
for gix = 1:numel(info.Groups)
    parts = strsplit(info.Groups(gix).Name, '/');
    grpName = parts{end};
    if startsWith(grpName, 'Path') || startsWith(grpName, 'DMD')
        pathNames{end+1} = grpName; %#ok<AGROW>
    end
end

if isempty(pathNames)
    ROIs = struct([]);
    return
end

pathNums = nan(size(pathNames));
for ix = 1:numel(pathNames)
    num = sscanf(pathNames{ix}, 'Path%d', 1);
    if isempty(num)
        num = sscanf(pathNames{ix}, 'DMD%d', 1);
    end
    if ~isempty(num)
        pathNums(ix) = num;
    end
end
[~, order] = sort(pathNums);
pathNames = pathNames(order);

ROIs = repmat(struct('dr', '', 'fn', '', 'roiData', {{}}, 'IM', [], ...
    'alignmentShiftRC', [], 'alignmentNCC', []), 1, numel(pathNames));
for DMDix = 1:numel(pathNames)
    pathRoot = ['/' pathNames{DMDix}];
    ROIs(DMDix).dr = readStringOrDefault(filename, [pathRoot '/dr'], '');
    ROIs(DMDix).fn = readStringOrDefault(filename, [pathRoot '/fn'], '');
    ROIs(DMDix).IM = readNumericOrDefault(filename, [pathRoot '/image'], []);
    ROIs(DMDix).alignmentShiftRC = readNumericOrDefault(filename, [pathRoot '/alignment_shift_rc'], []);
    ROIs(DMDix).alignmentNCC = readNumericOrDefault(filename, [pathRoot '/alignment_ncc'], []);

    nRois = double(readNumericOrDefault(filename, [pathRoot '/n_rois'], uint32(0)));
    roiData = cell(1, nRois);
    for rix = 1:nRois
        roiPath = sprintf('%s/roi_%03d', pathRoot, rix);

        typeRaw = readStringOrDefault(filename, [roiPath '/type'], '');
        S = struct();
        S.Type = denormalizeType(typeRaw);
        S.Label = readStringOrDefault(filename, [roiPath '/label'], '');
        S.mask = logical(readNumericOrDefault(filename, [roiPath '/mask'], uint8([])));

        position = readNumericOrDefault(filename, [roiPath '/position'], []);
        if ~isempty(position)
            S.Position = h5CoordsToMatlabRoi(double(position), coordsZeroIndexed);
        end

        center = readNumericOrDefault(filename, [roiPath '/center'], []);
        if ~isempty(center)
            S.Center = h5CoordsToMatlabRoi(double(center), coordsZeroIndexed);
        end

        semiAxes = readNumericOrDefault(filename, [roiPath '/semi_axes'], []);
        if ~isempty(semiAxes)
            S.SemiAxes = double(semiAxes(:).');
            S.AspectRatio = S.SemiAxes(1) ./ max(eps, S.SemiAxes(2));
            S.FixedAspectRatio = false;
        end

        rotationAngle = readNumericOrDefault(filename, [roiPath '/rotation_angle'], []);
        if ~isempty(rotationAngle)
            S.RotationAngle = double(rotationAngle(1));
        end

        radius = readNumericOrDefault(filename, [roiPath '/radius'], []);
        if ~isempty(radius)
            S.Radius = double(radius(1));
        end

        roiData{rix} = S;
    end
    ROIs(DMDix).roiData = roiData;
end
end

function value = readStringOrDefault(filename, path, defaultValue)
try
    raw = h5read(filename, path);
    value = char(uint16(raw(:)).');
catch
    value = defaultValue;
end
end

function value = readNumericOrDefault(filename, path, defaultValue)
try
    value = h5read(filename, path);
catch
    value = defaultValue;
end
end

function xy = h5CoordsToMatlabRoi(coords, zeroIndexed)
%H5COORDSTOMATLABROI Convert H5 coords to images.roi [x, y] (1-based).
coords = double(coords);
if zeroIndexed
    xy = [coords(:, 2) + 1, coords(:, 1) + 1];
else
    xy = coords;
end
if size(xy, 1) == 1
    xy = xy(1, :);
end
end


function typeOut = denormalizeType(typeIn)
switch lower(strtrim(typeIn))
    case 'polygon'
        typeOut = 'images.roi.polygon';
    case 'circle'
        typeOut = 'images.roi.circle';
    case 'ellipse'
        typeOut = 'images.roi.ellipse';
    otherwise
        typeOut = typeIn;
end
end
