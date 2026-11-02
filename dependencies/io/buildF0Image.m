function img = buildF0Image(E, selPix, numChannels, orderedChannels)
%BUILDF0IMAGE Per-trial mean background (F0) image from SILo trial results.
%
%   IMG = BUILDF0IMAGE(E, SELPIX, NUMCHANNELS, ORDEREDCHANNELS) scatters
%   E.F0im (selected pixels x analyzed channels, ordered channel order) into
%   a rows x cols x NUMCHANNELS single image in original channel order.
%   Pixels outside SELPIX and channels that were not analyzed are NaN.
%   Returns [] if E has no usable F0im.

img = [];
if isempty(E) || ~isfield(E, 'F0im') || isempty(E.F0im) || isempty(selPix)
    return
end
if size(E.F0im, 1) ~= nnz(selPix)
    warning('buildF0Image:SizeMismatch', ...
        'F0im has %d pixels but selPix has %d; skipping this trial.', ...
        size(E.F0im, 1), nnz(selPix));
    return
end
[imRows, imCols] = size(selPix);
nAnalyzed = size(E.F0im, 2);
img = nan(imRows*imCols, numChannels, 'single');
img(selPix(:), orderedChannels(1:nAnalyzed)) = E.F0im;
img = reshape(img, imRows, imCols, numChannels);
end
