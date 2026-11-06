function [Uout,S,Vout, BG] = nansvd(D, numComps, maxIter, nanFrac)
%performs sparse svd on data with nans by iteratively filling gaps with the best
%estimate of the requested rank

%for correct output, pixels should be in rows and timepoints in columns

if nargin<3
    maxIter = 20;
end
if nargin<4
    nanFrac = 0.3;
end

%select valid data to reconstruct
valid = ~isnan(D);
selCols = true(1,size(D,2));
selRows = true(1,size(D,1));
nC = 0; nR = 0; 
while ~(sum(selCols)==nC && sum(selRows)==nR)
    nC = sum(selCols);
    nR = sum(selRows);
    selRows = mean(~valid(:,selCols),2)<nanFrac;
    selCols = mean(~valid(selRows,:),1)<nanFrac;
end
D = D(selRows,selCols);
bg = smoothdata(D,2, 'movmean',length(selCols)/5);
D = D-bg;

nans = isnan(D);
estimate = repmat(mean(D,2, 'omitnan'),1, size(D,2));
for iter = 1:maxIter
   %fill nans with estimate
   D(nans) = estimate(nans);

   %perform svd
   [U,S,V] = svds(D,numComps);

   %update estimate
   estimate = U*S*V';

   %check convergence
   convCrit = max(((estimate(nans)-D(nans))./D(nans)).^2); %convergence criterion
    if convCrit<1e-6  %if converged, break
        break
    end
end

%svds returns min(numComps, min(size(D))) components, so fewer when D has
%fewer valid pixels or timepoints than numComps; pad so outputs always have
%numComps components (NaN in U and V, zero singular values in S)
nFound = size(U,2);
Uout = nan(length(selRows),numComps);
Uout(selRows,1:nFound) = U;
Vout = nan(length(selCols), numComps);
Vout(selCols,1:nFound) = V;
if nFound < numComps
    S(numComps,numComps) = 0; %grows S, filling the new entries with zeros
end
BG = nan(length(selRows),length(selCols));
BG(selRows,selCols) = bg;
end