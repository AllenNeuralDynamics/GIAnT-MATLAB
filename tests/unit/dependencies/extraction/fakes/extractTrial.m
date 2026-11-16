function future = extractTrial(Y_obs, Finv, sources, selPix, params, alignData) %#ok<INUSD>
%EXTRACTTRIAL Test double: saves the samples processAllTrials_Async hands to
%source extraction (pixels x frames x channels) to the file named by
%getappdata(0, 'FakeExtractTrialFile'), and returns a future with empty outputs.
fn = getappdata(0, 'FakeExtractTrialFile');
if ~isempty(fn)
    save(fn, 'Y_obs', 'Finv', 'selPix');
end
future = FakeExtractFuture();
end
