function [motion, R] = xcorr2_nans_weighted(frame, freshness, template, shiftsCenter, dShift)
%Kaspar Podgorski 2023
%perform a somewhat-efficient local normalized crosscorrelation for images with
%nans

%template: the template
%frame: the frame to be aligned; this has more NaNs

%shiftsCenter: the center offset around which to perform a local search
%dShift: the maximum shift (scalar, in pixels) to consider on each axis around shiftsCenter
%freshness: the weight of each frame pixel. For each shift, the correlation
%is the freshness-weighted Pearson correlation over the pixels valid in both
%images: weighted means, variances and covariance, so that a pixel's influence
%is set by its weight alone. (Mixing unweighted means and template variance
%with weighted sums biased the peak when the set of compared pixels changed
%with the shift, e.g. on thin raster ribbons.)

sz= size(template);
dShift = round(dShift(1)); %sanity check

template = circshift(template, shiftsCenter);

fIdx = find(~isnan(frame)); %valid pixels of the fixed image, column-major
[Fr0,Fc0] = ind2sub(sz, fIdx);
F0 = frame(fIdx); Ff0 = freshness(fIdx);

%correlation is sum(A.*B)./(sqrt(ssA)*sqrt(ssB)); ssB is constant though
%For each shift, the valid frame pixels whose shifted template pixel is inside
%the image and not NaN, in column-major order; their template pixels are at
%linear index fIdx + dr + dc*sz(1)
shifts =  -dShift:dShift;
inC = cell(size(shifts));
for dcix = 1:length(shifts)
    inC{dcix} = Fc0+shifts(dcix)>=1 & Fc0+shifts(dcix)<=sz(2);
end
C = nan(length(shifts), length(shifts));
for drix = 1:length(shifts)
    inR = Fr0+shifts(drix)>=1 & Fr0+shifts(drix)<=sz(1);
    for dcix = 1:length(shifts)
        valid = inR & inC{dcix};
        T = template(fIdx(valid) + (shifts(drix) + shifts(dcix)*sz(1)));
        sel = ~isnan(T); T = T(sel);
        F = F0(valid); F = F(sel);
        Ff = Ff0(valid); Ff = Ff(sel);

        sw = sum(Ff);
        dF = F - sum(Ff.*F)/sw; dT = T - sum(Ff.*T)/sw; %weighted means
        sFT = sum(Ff.*dF.*dT);
        sT = sum(Ff.*dT.^2);
        sF = sum(Ff.*dF.^2);
        C(drix,dcix) = sFT./sqrt(sT.*sF);
    end
end

%find maximum of correlation map
[maxval,I] = max(C(:));
[rr,cc] = ind2sub(size(C),I);

R = maxval; %correlation coefficient

if rr>1 && rr<length(shifts) && cc>1 && cc<length(shifts)
    %perform superresolution upsampling
    ratioR = min(1e6,(C(rr,cc) - C(rr-1,cc))/(C(rr,cc) - C(rr+1,cc)));
    dR = (1-ratioR)/(1+ratioR)/2;
    
    ratioC =min(1e6, (C(rr,cc) - C(rr,cc-1))/(C(rr,cc) - C(rr,cc+1)));
    dC = (1-ratioC)/(1+ratioC)/2;

    motion = -shiftsCenter' + [shifts(rr)-dR shifts(cc)-dC];
else %the optimum is at an edge of search range; no superresolution
    motion = -shiftsCenter' + [shifts(rr) shifts(cc)];
end

motion = -motion;

if any(isnan(motion))
    error('xcorr2_nans_weighted:NaNMotion', 'Computed motion contains NaN.');
end
