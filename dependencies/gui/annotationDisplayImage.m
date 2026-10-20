function IMdisp = annotationDisplayImage(IM)
%ANNOTATIONDISPLAYIMAGE Image shown by drawROIs, scaled to [0 1].
%
% IM is H x W (grayscale) or H x W x 2 (shown as ch1 red, ch2 cyan). Each
% channel is divided by its 99.9th percentile over imaged pixels (finite and
% > 0). Sparse SLAP2 images are mostly unimaged background, so a percentile
% over all pixels can be 0, which would divide by zero and display black.

switch size(IM,3)
    case 1
        IMdisp = scaleChannel(IM);
    case 2
        IM1 = scaleChannel(IM(:,:,1));
        IM2 = scaleChannel(IM(:,:,2));
        IMdisp = cat(3, IM1, IM2, IM2);
    otherwise
        error('Unexpected number of channels (dim3) for drawROIs GUI');
end
IMdisp = min(max(IMdisp, 0), 1);
IMdisp(isnan(IMdisp)) = 0;
end


function x = scaleChannel(x)
x = double(x);
vals = x(isfinite(x) & x > 0);
if isempty(vals)
    return %nothing imaged in this channel; leave it black
end
x = x ./ prctile(vals, 99.9);
end
