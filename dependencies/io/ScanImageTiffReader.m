classdef ScanImageTiffReader < handle
    %SCANIMAGETIFFREADER Pure MATLAB reader for ScanImage TIFF files.
    %   ScanImageTiffReader that does not require
    %   MEX files. Reads uncompressed strip-based TIFF and BigTIFF files
    %   written by ScanImage.
    %
    %   Usage:
    %       reader = ScanImageTiffReader('path/to/file.tif');
    %       vol    = reader.data();
    %       descs  = reader.descriptions();
    %       meta   = reader.metadata();
    %       reader.close();

    properties (Access = private)
        filename_
        fid_        = -1
        machfmt_             % 'ieee-le' or 'ieee-be'
        needSwap_   = false  % true when file byte order differs from native
        isBigTiff_  = false

        firstIfdOffset_

        % Per-IFD index (cell arrays, one element per IFD / frame)
        ifdTags_             % cell of struct arrays:  .id .type .count .dataOffset .dataNbytes
        ifdStrips_           % cell of struct arrays:  .offset .nbytes

        % Image geometry (determined from first IFD)
        w_          = 0
        h_          = 0
        nchan_      = 1
        nplanes_    = 0
        matlabType_          % e.g. 'uint16', 'single'
        ndType_              % integer matching the C enum nd_type
        bytesPerSample_ = 0

        endOfImageBlock_ = 0
        fileSize_        = 0
    end

    % -----------------------------------------------------------------
    methods (Static)
        function out = apiVersion()
            out = 'ScanImageTiffReader (pure MATLAB)';
        end

        function buf = readByteRange_(filename, offset, nbytes)
            %READBYTERANGE_ Read nbytes starting at offset, with a private file handle.
            fid = fopen(filename, 'r');
            assert(fid > 0, 'ScanImageTiffReader:OpenFailed', 'Could not open file %s', filename);
            closeFile = onCleanup(@() fclose(fid));
            fseek(fid, offset, 'bof');
            buf = fread(fid, nbytes, '*uint8');
            assert(numel(buf) == nbytes, 'ScanImageTiffReader:ShortRead', ...
                'Read %d of %d bytes at offset %d of %s', numel(buf), nbytes, offset, filename);
        end

    end

    % -----------------------------------------------------------------
    methods
        function obj = ScanImageTiffReader(filename)
            if nargin > 0
                obj.open(filename);
            end
        end

        function delete(obj)
            obj.close();
        end

        function obj = open(obj, filename)
            if obj.fid_ > 0
                obj.close();
            end

            assert(exist(filename, 'file') == 2, ...
                'ScanImageTiffReader:FileNotFound', ...
                'File %s not found on disk.', filename);

            listing = dir(filename);
            obj.filename_ = fullfile(listing.folder, listing.name);

            obj.fid_ = fopen(obj.filename_, 'r');
            assert(obj.fid_ > 0, ...
                'ScanImageTiffReader:OpenFailed', ...
                'Could not open file %s', obj.filename_);

            fseek(obj.fid_, 0, 'eof');
            obj.fileSize_ = ftell(obj.fid_);

            obj.identify_();
            obj.buildIndex_();
        end

        function obj = close(obj)
            if obj.fid_ > 0
                fclose(obj.fid_);
                obj.fid_ = -1;
            end
        end

        function tf = isOpen(obj)
            tf = obj.fid_ > 0;
        end

        % ----- data --------------------------------------------------
        function stack = data(obj)
            obj.ensureOpen_();

            nTotalStrips = 0;
            for k = 1:obj.nplanes_
                nTotalStrips = nTotalStrips + numel(obj.ifdStrips_{k});
            end

            allOffsets = zeros(nTotalStrips, 1);
            allNbytes  = zeros(nTotalStrips, 1);
            idx = 0;
            for k = 1:obj.nplanes_
                ss = obj.ifdStrips_{k};
                for j = 1:numel(ss)
                    idx = idx + 1;
                    allOffsets(idx) = ss(j).offset;
                    allNbytes(idx)  = ss(j).nbytes;
                end
            end

            totalBytes   = sum(allNbytes);
            totalSamples = totalBytes / obj.bytesPerSample_;
            precision    = ['*' obj.matlabType_];

            % Fast path: single read when all strips are contiguous
            contiguous = true;
            expected   = allOffsets(1);
            for k = 1:nTotalStrips
                if allOffsets(k) ~= expected
                    contiguous = false;
                    break;
                end
                expected = expected + allNbytes(k);
            end

            if contiguous
                fseek(obj.fid_, allOffsets(1), 'bof');
                stack = fread(obj.fid_, totalSamples, precision, 0, obj.machfmt_);
            else
                stack    = zeros(totalSamples, 1, obj.matlabType_);
                samplePos = 1;
                for k = 1:nTotalStrips
                    nSamp = allNbytes(k) / obj.bytesPerSample_;
                    fseek(obj.fid_, allOffsets(k), 'bof');
                    stack(samplePos:samplePos + nSamp - 1) = ...
                        fread(obj.fid_, nSamp, precision, 0, obj.machfmt_);
                    samplePos = samplePos + nSamp;
                end
            end

            dims = [obj.nchan_, obj.w_, obj.h_, obj.nplanes_];
            dims(dims == 1) = [];
            if numel(dims) > 1
                stack = reshape(stack, dims);
            end
        end

        % ----- dataPlanes --------------------------------------------
        function stack = dataPlanes(obj, planeIdx, nParallel)
            %DATAPLANES Read only the given planes (IFDs), in the order given.
            %   stack = reader.dataPlanes(planeIdx) returns the same layout as
            %   data() with nplanes = numel(planeIdx), reading each plane's
            %   strips directly, so memory scales with the planes requested
            %   rather than the file. Nearby planes are read together in
            %   sequential reads of up to 64 MB.
            %
            %   stack = reader.dataPlanes(planeIdx, nParallel) issues up to
            %   nParallel of those reads concurrently (thread-based
            %   backgroundPool). This is much faster on storage where each
            %   read waits on the network (object-store mounts, NFS, SMB).
            %   Falls back to sequential reads if no pool is available.
            if nargin < 3 || isempty(nParallel)
                nParallel = 1;
            end
            obj.ensureOpen_();
            planeIdx = planeIdx(:)';
            assert(all(planeIdx >= 1 & planeIdx <= obj.nplanes_ & planeIdx == round(planeIdx)), ...
                'ScanImageTiffReader:InvalidPlane', 'Plane indices must be integers in 1..%d.', obj.nplanes_);

            samplesPerPlane = obj.nchan_ * obj.w_ * obj.h_;
            stack = zeros(samplesPerPlane, numel(planeIdx), obj.matlabType_);
            [runs, scattered] = obj.planReads_(planeIdx, samplesPerPlane);

            % planes whose strips are scattered within the plane: read strip by strip
            precision = ['*' obj.matlabType_];
            for p = scattered
                ss = obj.ifdStrips_{planeIdx(p)};
                pos = 1;
                for j = 1:numel(ss)
                    nSamp = ss(j).nbytes / obj.bytesPerSample_;
                    fseek(obj.fid_, ss(j).offset, 'bof');
                    stack(pos:pos + nSamp - 1, p) = fread(obj.fid_, nSamp, precision, 0, obj.machfmt_);
                    pos = pos + nSamp;
                end
            end

            if nParallel > 1 && numel(runs) > 1
                try
                    stack = obj.readRunsParallel_(runs, stack, nParallel);
                    runs = [];
                catch err
                    warning('ScanImageTiffReader:ParallelReadFailed', ...
                        'Parallel read failed (%s); reading sequentially.', err.message);
                end
            end
            for r = 1:numel(runs)
                fseek(obj.fid_, runs(r).start, 'bof');
                buf = fread(obj.fid_, runs(r).stop - runs(r).start, '*uint8');
                stack = obj.unpackRun_(runs(r), buf, stack);
            end

            dims = [obj.nchan_, obj.w_, obj.h_, numel(planeIdx)];
            dims(dims == 1) = [];
            if numel(dims) > 1
                stack = reshape(stack, dims);
            end
        end

        % ----- compression -------------------------------------------
        function c = compression(obj)
            %COMPRESSION TIFF Compression tag of the first IFD (1 = none).
            %   This reader only decodes uncompressed data.
            obj.ensureOpen_();
            tags = obj.ifdTags_{1};
            if any([tags.id] == 259)
                c = double(obj.readTagValue_(tags, 259));
            else
                c = 1;
            end
        end

        % ----- descriptions ------------------------------------------
        function desc = descriptions(obj)
            obj.ensureOpen_();
            desc = cell(obj.nplanes_, 1);
            for k = 1:obj.nplanes_
                tags = obj.ifdTags_{k};
                ti = find([tags.id] == 270, 1);   % IMAGEDESCRIPTION
                if ~isempty(ti)
                    tag = tags(ti);
                    fseek(obj.fid_, tag.dataOffset, 'bof');
                    raw = fread(obj.fid_, tag.count, '*char')';
                    ni = find(raw == 0, 1);
                    if ~isempty(ni)
                        raw = raw(1:ni - 1);
                    end
                    desc{k} = strtrim(raw);
                else
                    desc{k} = '';
                end
            end
        end

        % ----- metadata ----------------------------------------------
        function meta = metadata(obj)
            obj.ensureOpen_();
            fmt = obj.detectMetadataFormat_();
            switch fmt
                case 1,     meta = obj.readMetadata2016_();
                case 2,     meta = obj.readMetadataOld_();
                otherwise,  meta = '';
            end
        end

        % ----- shape -------------------------------------------------
        function s = shape(obj)
            obj.ensureOpen_();
            dims = [obj.nchan_, obj.w_, obj.h_, obj.nplanes_];
            dims(dims == 1) = [];

            s.ndim    = numel(dims);
            s.dims    = dims;
            s.type    = obj.ndType_;
            s.strides = zeros(1, s.ndim + 1);
            s.strides(1) = 1;
            for k = 1:s.ndim
                s.strides(k + 1) = dims(k) * s.strides(k);
            end
        end
    end

    % =================================================================
    %  Private helpers
    % =================================================================
    methods (Access = private)

        % ----- planReads_ -------------------------------------------
        function [runs, scattered] = planReads_(obj, planeIdx, samplesPerPlane)
            % Group the requested planes into sequential reads ("runs"):
            % nearby planes stored as contiguous blocks are read together,
            % reading through gaps up to maxGapBytes, each run at most
            % maxRunBytes. Planes with scattered strips are returned
            % separately (positions in planeIdx).
            maxGapBytes = 8e6;
            maxRunBytes = 64e6;
            nReq = numel(planeIdx);
            planeStart = zeros(1, nReq);
            planeEnd   = zeros(1, nReq);
            isBlock    = false(1, nReq);
            for p = 1:nReq
                ss = obj.ifdStrips_{planeIdx(p)};
                offsets = [ss.offset];
                nbytes  = [ss.nbytes];
                isBlock(p) = all(offsets(2:end) == offsets(1:end-1) + nbytes(1:end-1)) ...
                    && sum(nbytes) == samplesPerPlane * obj.bytesPerSample_;
                planeStart(p) = offsets(1);
                planeEnd(p)   = offsets(1) + sum(nbytes);
            end
            scattered = find(~isBlock);
            blocks = find(isBlock);
            [~, o] = sort(planeStart(blocks));
            blocks = blocks(o);

            runs = struct('start', {}, 'stop', {}, 'planes', {}, 'planeStart', {}, 'planeEnd', {});
            k = 1;
            while k <= numel(blocks)
                runStart = planeStart(blocks(k));
                runEnd   = planeEnd(blocks(k));
                last = k;
                while last < numel(blocks)
                    q = blocks(last + 1);
                    if planeStart(q) - runEnd > maxGapBytes || planeEnd(q) - runStart > maxRunBytes
                        break
                    end
                    runEnd = max(runEnd, planeEnd(q));
                    last = last + 1;
                end
                members = blocks(k:last);
                runs(end+1) = struct('start', runStart, 'stop', runEnd, 'planes', members, ...
                    'planeStart', planeStart(members), 'planeEnd', planeEnd(members)); %#ok<AGROW>
                k = last + 1;
            end
        end

        % ----- unpackRun_ -------------------------------------------
        function stack = unpackRun_(obj, run, buf, stack)
            for m = 1:numel(run.planes)
                a = run.planeStart(m) - run.start + 1; %scalar range indexing copies a contiguous block
                b = run.planeEnd(m) - run.start;
                vals = typecast(buf(a:b), obj.matlabType_);
                if obj.needSwap_
                    vals = swapbytes(vals);
                end
                stack(:, run.planes(m)) = vals;
            end
        end

        % ----- readRunsParallel_ ------------------------------------
        function stack = readRunsParallel_(obj, runs, stack, nParallel)
            % Each run is read by a backgroundPool thread with its own file
            % handle; at most nParallel reads are in flight, and runs are
            % unpacked as they complete so buffers do not accumulate.
            pool = backgroundPool;
            nParallel = min(nParallel, pool.NumWorkers);
            if nParallel < 2
                error('ScanImageTiffReader:NoParallelPool', 'backgroundPool has %d worker(s).', pool.NumWorkers);
            end
            fn = obj.filename_;
            futures = parallel.FevalFuture.empty;
            owner = zeros(1, 0);
            next = 1;
            nDone = 0;
            cancelAll = onCleanup(@() cancel(futures));
            while nDone < numel(runs)
                while next <= numel(runs) && sum(strcmp({futures.State}, 'queued') | strcmp({futures.State}, 'running')) < nParallel
                    futures(end+1) = parfeval(pool, @ScanImageTiffReader.readByteRange_, 1, ...
                        fn, runs(next).start, runs(next).stop - runs(next).start); %#ok<AGROW>
                    owner(end+1) = next; %#ok<AGROW>
                    next = next + 1;
                end
                [idx, buf] = fetchNext(futures);
                stack = obj.unpackRun_(runs(owner(idx)), buf, stack);
                nDone = nDone + 1;
            end
        end


        function ensureOpen_(obj)
            if obj.fid_ < 0
                error('ScanImageTiffReader:NotOpen', ...
                    'File is not open for reading.');
            end
        end

        % ----- identify_ --------------------------------------------
        function identify_(obj)
            fseek(obj.fid_, 0, 'bof');
            bom = fread(obj.fid_, 2, '*uint8');

            if isequal(bom, uint8([73; 73]))       % 'II'
                obj.machfmt_ = 'ieee-le';
            elseif isequal(bom, uint8([77; 77]))   % 'MM'
                obj.machfmt_ = 'ieee-be';
            else
                error('ScanImageTiffReader:InvalidFile', ...
                    'Invalid TIFF: unrecognized byte order mark.');
            end

            [~, ~, ne] = computer;
            fileIsLE   = strcmp(obj.machfmt_, 'ieee-le');
            nativeIsLE = (ne == 'L');
            obj.needSwap_ = (fileIsLE ~= nativeIsLE);

            subtype = fread(obj.fid_, 1, 'uint16', 0, obj.machfmt_);
            switch subtype
                case 42
                    obj.isBigTiff_ = false;
                case 43
                    obj.isBigTiff_ = true;
                otherwise
                    error('ScanImageTiffReader:InvalidFile', ...
                        'Unrecognized TIFF sub-type %d.', subtype);
            end

            if obj.isBigTiff_
                fseek(obj.fid_, 8, 'bof');
                obj.firstIfdOffset_ = fread(obj.fid_, 1, 'uint64', 0, obj.machfmt_);
            else
                obj.firstIfdOffset_ = fread(obj.fid_, 1, 'uint32', 0, obj.machfmt_);
            end
        end

        % ----- buildIndex_ ------------------------------------------
        function buildIndex_(obj)
            obj.ifdTags_   = {};
            obj.ifdStrips_ = {};

            if obj.isBigTiff_
                entrySz       = 20;
                nentriesFmt   = 'uint64';
                offsetFmt     = 'uint64';
                dataFieldOff  = 12;     % offset within entry to the data/pointer field
                dataFieldSz   = 8;
            else
                entrySz       = 12;
                nentriesFmt   = 'uint16';
                offsetFmt     = 'uint32';
                dataFieldOff  = 8;
                dataFieldSz   = 4;
            end

            offset   = obj.firstIfdOffset_;
            planeIdx = 0;

            while offset ~= 0
                planeIdx = planeIdx + 1;

                fseek(obj.fid_, offset, 'bof');
                nentries   = fread(obj.fid_, 1, nentriesFmt, 0, obj.machfmt_);
                entryStart = ftell(obj.fid_);

                raw = fread(obj.fid_, nentries * entrySz, '*uint8');

                nextOffset = fread(obj.fid_, 1, offsetFmt, 0, obj.machfmt_);

                tags = obj.parseEntries_(raw, nentries, entrySz, ...
                    dataFieldOff, dataFieldSz, entryStart);

                strips = obj.readStripInfo_(tags);

                if planeIdx == 1
                    obj.w_ = obj.readTagValue_(tags, 256);    % IMAGEWIDTH
                    obj.h_ = obj.readTagValue_(tags, 257);    % IMAGELENGTH

                    if any([tags.id] == 277)
                        obj.nchan_ = obj.readTagValue_(tags, 277);  % SAMPLESPERPIXEL
                    else
                        obj.nchan_ = 1;
                    end

                    bps = obj.readTagValue_(tags, 258);       % BITSPERSAMPLE
                    obj.bytesPerSample_ = bps / 8;

                    sf = 1;  % default: unsigned
                    if any([tags.id] == 339)
                        sf = obj.readTagValue_(tags, 339);    % SAMPLEFORMAT
                    end

                    [obj.matlabType_, obj.ndType_] = ...
                        ScanImageTiffReader.resolvePixelType_(sf, bps);
                end

                obj.ifdTags_{planeIdx}   = tags;
                obj.ifdStrips_{planeIdx} = strips;

                offset = nextOffset;
            end

            obj.nplanes_ = planeIdx;

            maxEnd = 0;
            for k = 1:obj.nplanes_
                ss = obj.ifdStrips_{k};
                for j = 1:numel(ss)
                    maxEnd = max(maxEnd, ss(j).offset + ss(j).nbytes);
                end
            end
            obj.endOfImageBlock_ = maxEnd;
        end

        % ----- parseEntries_ ----------------------------------------
        function tags = parseEntries_(obj, raw, nentries, entrySz, ...
                dataFieldOff, dataFieldSz, entryStart)
            sw = obj.needSwap_;
            bt = obj.isBigTiff_;

            tags = repmat(struct('id',0,'type',0,'count',0, ...
                'dataOffset',0,'dataNbytes',0), 1, nentries);

            for i = 1:nentries
                b = (i - 1) * entrySz;

                tagId  = typecast(raw(b+1  : b+2),  'uint16');
                typeId = typecast(raw(b+3  : b+4),  'uint16');
                if bt
                    cnt = typecast(raw(b+5  : b+12), 'uint64');
                    df  = typecast(raw(b+13 : b+20), 'uint64');
                else
                    cnt = typecast(raw(b+5  : b+8),  'uint32');
                    df  = typecast(raw(b+9  : b+12), 'uint32');
                end

                if sw
                    tagId  = swapbytes(tagId);
                    typeId = swapbytes(typeId);
                    cnt    = swapbytes(cnt);
                    df     = swapbytes(df);
                end

                tsz       = ScanImageTiffReader.tiffTypeSize_(double(typeId));
                dataNb    = double(cnt) * tsz;

                if dataNb <= dataFieldSz
                    doff = entryStart + b + dataFieldOff;
                else
                    doff = double(df);
                end

                tags(i).id         = double(tagId);
                tags(i).type       = double(typeId);
                tags(i).count      = double(cnt);
                tags(i).dataOffset = doff;
                tags(i).dataNbytes = dataNb;
            end
        end

        % ----- readTagValue_ ----------------------------------------
        function val = readTagValue_(obj, tags, tagId)
            ti = find([tags.id] == tagId, 1);
            assert(~isempty(ti), ...
                'ScanImageTiffReader:TagNotFound', ...
                'Required tag %d not found.', tagId);
            tag = tags(ti);
            fseek(obj.fid_, tag.dataOffset, 'bof');
            mtype = ScanImageTiffReader.tiffTypeToMatlab_(tag.type);
            val = fread(obj.fid_, 1, mtype, 0, obj.machfmt_);
        end

        % ----- readStripInfo_ ---------------------------------------
        function strips = readStripInfo_(obj, tags)
            oIdx = find([tags.id] == 273, 1);   % STRIPOFFSETS
            cIdx = find([tags.id] == 279, 1);   % STRIPBYTECOUNTS

            if isempty(oIdx) || isempty(cIdx)
                error('ScanImageTiffReader:NoStrips', ...
                    'No strip offsets / byte counts found. Only strip-based TIFFs are supported.');
            end

            oTag = tags(oIdx);
            cTag = tags(cIdx);

            oType = ScanImageTiffReader.tiffTypeToMatlab_(oTag.type);
            cType = ScanImageTiffReader.tiffTypeToMatlab_(cTag.type);

            fseek(obj.fid_, oTag.dataOffset, 'bof');
            offsets = fread(obj.fid_, oTag.count, oType, 0, obj.machfmt_)';

            fseek(obj.fid_, cTag.dataOffset, 'bof');
            nbytes = fread(obj.fid_, cTag.count, cType, 0, obj.machfmt_)';

            strips = struct('offset', num2cell(offsets), ...
                            'nbytes', num2cell(nbytes));
        end

        % ----- metadata format detection ----------------------------
        function fmt = detectMetadataFormat_(obj)
            % Returns 0 = none, 1 = SI 2016+, 2 = legacy (appended block)
            if obj.fileSize_ < 32
                fmt = 0;
                return;
            end

            fseek(obj.fid_, 16, 'bof');
            magic = fread(obj.fid_, 1, 'uint32', 0, obj.machfmt_);

            if magic == hex2dec('07030301')
                fmt = 1;
                return;
            end

            if obj.endOfImageBlock_ ~= obj.fileSize_
                fmt = 2;
                return;
            end

            fmt = 0;
        end

        % ----- readMetadata2016_ ------------------------------------
        function meta = readMetadata2016_(obj)
            fseek(obj.fid_, 24, 'bof');   % skip magic + versionId
            headerBytes   = fread(obj.fid_, 1, 'uint32', 0, obj.machfmt_);
            roigroupBytes = fread(obj.fid_, 1, 'uint32', 0, obj.machfmt_);

            nbytes = headerBytes + roigroupBytes;
            if nbytes == 0
                meta = '';
                return;
            end

            fseek(obj.fid_, 32, 'bof');
            buf = fread(obj.fid_, nbytes, '*char')';

            if roigroupBytes == 0
                meta = deblank(buf);
                return;
            end

            if buf(1) == '{'
                a = ScanImageTiffReader.nullTerminate_(buf(1:headerBytes));
                b = ScanImageTiffReader.nullTerminate_(buf(headerBytes+1:end));
                a = ScanImageTiffReader.removeEnclosingBraces_(a);
                b = ScanImageTiffReader.removeEnclosingBraces_(b);
                meta = ['{' a ',' b '}'];
            else
                buf(headerBytes) = char(10);   % replace null separator with newline
                meta = buf;
            end
        end

        % ----- readMetadataOld_ -------------------------------------
        function meta = readMetadataOld_(obj)
            fseek(obj.fid_, obj.fileSize_ - 4, 'bof');
            nbytes = fread(obj.fid_, 1, 'uint32', 0, obj.machfmt_);

            if nbytes == 0 || nbytes > obj.fileSize_
                meta = '';
                return;
            end

            fseek(obj.fid_, obj.fileSize_ - 4 - nbytes, 'bof');
            meta = fread(obj.fid_, nbytes, '*char')';
        end
    end

    % =================================================================
    %  Static private helpers (no file access needed)
    % =================================================================
    methods (Static, Access = private)

        function sz = tiffTypeSize_(typeId)
            persistent sizes
            if isempty(sizes)
                sizes = zeros(1, 18);
                %       BYTE ASCII SHORT LONG  RATIONAL SBYTE UNDEFINED SSHORT SLONG SRATIONAL FLOAT DOUBLE IFD4       LONG8 SLONG8 IFD8
                sizes = [1    1     2     4     8        1     0         2      4     8         4     8      4    0 0    8     8      8   ];
            end
            if typeId >= 1 && typeId <= numel(sizes)
                sz = sizes(typeId);
            else
                error('ScanImageTiffReader:BadType', ...
                    'Unrecognized TIFF element type %d.', typeId);
            end
        end

        function mtype = tiffTypeToMatlab_(typeId)
            switch typeId
                case {1,2,7},   mtype = 'uint8';    % BYTE / ASCII / UNDEFINED
                case 3,         mtype = 'uint16';   % SHORT
                case {4,13},    mtype = 'uint32';   % LONG / IFD4
                case 5,         mtype = 'uint32';   % RATIONAL (read as raw uint32 pairs)
                case 6,         mtype = 'int8';     % SBYTE
                case 8,         mtype = 'int16';    % SSHORT
                case {9,10},    mtype = 'int32';    % SLONG / SRATIONAL
                case 11,        mtype = 'single';   % FLOAT
                case 12,        mtype = 'double';   % DOUBLE
                case {16,18},   mtype = 'uint64';   % LONG8 / IFD8
                case 17,        mtype = 'int64';    % SLONG8
                otherwise
                    error('ScanImageTiffReader:BadType', ...
                        'Unrecognized TIFF element type %d.', typeId);
            end
        end

        function [mtype, ndtype] = resolvePixelType_(sampleFormat, bitsPerSample)
            switch sampleFormat
                case 3   % FLOATING
                    switch bitsPerSample
                        case 32, mtype = 'single'; ndtype = 8;
                        case 64, mtype = 'double'; ndtype = 9;
                        otherwise, error('ScanImageTiffReader:BadPixel', ...
                                'Unsupported float bits per sample: %d', bitsPerSample);
                    end
                case 2   % SIGNED INTEGER
                    switch bitsPerSample
                        case 8,  mtype = 'int8';  ndtype = 4;
                        case 16, mtype = 'int16'; ndtype = 5;
                        case 32, mtype = 'int32'; ndtype = 6;
                        case 64, mtype = 'int64'; ndtype = 7;
                        otherwise, error('ScanImageTiffReader:BadPixel', ...
                                'Unsupported signed integer bits per sample: %d', bitsPerSample);
                    end
                otherwise  % UNSIGNED (1) or UNDEFINED (4)
                    switch bitsPerSample
                        case 8,  mtype = 'uint8';  ndtype = 0;
                        case 16, mtype = 'uint16'; ndtype = 1;
                        case 32, mtype = 'uint32'; ndtype = 2;
                        case 64, mtype = 'uint64'; ndtype = 3;
                        otherwise, error('ScanImageTiffReader:BadPixel', ...
                                'Unsupported unsigned integer bits per sample: %d', bitsPerSample);
                    end
            end
        end

        function s = nullTerminate_(s)
            ni = find(s == 0, 1);
            if ~isempty(ni)
                s = s(1:ni - 1);
            end
        end

        function s = removeEnclosingBraces_(s)
            bi = find(s == '{', 1, 'first');
            ei = find(s == '}', 1, 'last');
            if ~isempty(bi) && ~isempty(ei) && ei > bi + 1
                s = s(bi + 1 : ei - 1);
            end
        end
    end
end
