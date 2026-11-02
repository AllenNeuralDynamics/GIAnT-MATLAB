function [data, meta, desc] = ScanImageTiffWrapper(~)
%SCANIMAGETIFFWRAPPER Test double returning a small 2-channel, 3-plane
%   reference stack with SLAP2-style JSON frame descriptions.
channels = [1 2];
Zs = [0 1 2];

data = zeros(8, 8, numel(channels)*numel(Zs), 'int16');
meta = '';
desc = {};
for z = Zs
    for c = channels
        desc{end+1} = jsonencode(struct('channel', c, 'z', z, ...
            'dmdPixel2SampleTransform', eye(3))); %#ok<AGROW>
    end
end
end
