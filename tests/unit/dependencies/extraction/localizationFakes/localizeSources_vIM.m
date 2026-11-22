function [skIm, P] = localizeSources_vIM(IM, vIM, params, ~)
% Fake for TestLoadAndProcessTrialAsync: appends its inputs to
% getappdata(groot, 'localizeCalls') and fails the first
% getappdata(groot, 'localizeFailures') calls.
calls = getappdata(groot, 'localizeCalls');
calls{end+1} = struct('IM', IM, 'vIM', vIM, 'alignHz', params.alignHz);
setappdata(groot, 'localizeCalls', calls);
if numel(calls) <= getappdata(groot, 'localizeFailures')
    error('localizeSources_vIM:fake', 'Localization failed');
end
skIm = mean(IM, 3, 'omitnan');
P = struct('row', 2, 'col', 3, 'val', 1);
end
