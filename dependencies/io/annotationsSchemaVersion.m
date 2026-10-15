function v = annotationsSchemaVersion()
%ANNOTATIONSSCHEMAVERSION Schema version written to annotations H5 files.
%
% Increment when the layout or meaning of annotations.h5 /
% postacq_annotations.h5 changes in a way readers must know about, and
% update loadAnnotationsH5 / importPostAcqAnnotations to handle older files.
%   1 - unversioned files (no /schema_version dataset)
%   2 - adds /schema_version and optional per-path image, alignment_shift_rc
%       and alignment_ncc; post-acquisition annotations require image
v = 2;
end
