classdef Fast_BigTiff_Write < handle
    %FAST_BIGTIFF_WRITE Test double for the BigTIFF writer.
    %   Appends each page as raw single-precision values, so tests can count
    %   pages from the file size (4 bytes per pixel).

    properties (Access = private)
        fid
    end

    methods
        function obj = Fast_BigTiff_Write(filename, ~, ~)
            obj.fid = fopen(filename, 'w');
        end

        function WriteIMG(obj, img)
            fwrite(obj.fid, img, 'single');
        end

        function close(obj)
            fclose(obj.fid);
        end
    end
end
