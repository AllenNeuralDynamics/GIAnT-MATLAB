function d = subproblemSaveDir(params)
%SUBPROBLEMSAVEDIR the folder in which extraction saves its subproblems.
%   d = SUBPROBLEMSAVEDIR(params) is params.saveSubproblemsDir as text, or ''
%   when saving is off: the field missing or empty, blank text, false or 0
%   (as JSON parameters may give it). Any other value that is not text is an
%   error. See extractTrial.
d = '';
if ~isfield(params, 'saveSubproblemsDir')
    return
end
v = params.saveSubproblemsDir;
if isempty(v) || ((islogical(v) || isnumeric(v)) && isscalar(v) && ~v)
    return
end
if ~(ischar(v) || (isstring(v) && isscalar(v)))
    error('extractTrial:saveSubproblemsDir', ...
        'params.saveSubproblemsDir must be a folder name, or empty to save nothing.');
end
d = strtrim(char(v));
end
