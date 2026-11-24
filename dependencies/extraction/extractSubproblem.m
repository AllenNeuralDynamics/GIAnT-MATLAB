function varargout = extractSubproblem(Y_p, F_inv_p, sources_p, selPix_p, params)
%EXTRACTSUBPROBLEM Source extraction of one subproblem, as extractTrial solves it.
%   [H, S, F0, LS, SNR, errFinal, B] = EXTRACTSUBPROBLEM(Y_p, F_inv_p, sources_p, selPix_p, params)
%   runs extractTrial's subproblem solver on one subproblem's inputs, for
%   example those that extractTrial saves with params.saveSubproblemsDir set
%   (sub_###.mat: load it and pass its Y_p, F_inv_p, sources_p, selPix_p and
%   params). The inputs are in the photon-normalized units extractTrial
%   fits in (the trial's data divided by params.photonScale), and so are the
%   outputs:
%     H    footprints, pixels x sources (each sums to 1)
%     S    activity events, sources x time x channels; the denoised activity
%          is convn(S(:,:,1), params.k, 'same') (channel 2: params.k2)
%     F0   footprint-weighted baseline of each source, sources x time x channels
%     LS   least-squares activity, sources x time x channels
%     SNR  each source's SNR
%     B    baseline, pixels x time x channels
%   errFinal is NaN (it is set only for fits with ground truth).
%
%   This is a developer entry point (subproblemRefitTest uses it); SILo
%   does not call it. See also extractTrial, subproblemRefitTest.
solver = extractTrial('subproblemSolver');
[varargout{1:max(1, nargout)}] = solver(double(Y_p), double(F_inv_p), sources_p, selPix_p, params);
end
