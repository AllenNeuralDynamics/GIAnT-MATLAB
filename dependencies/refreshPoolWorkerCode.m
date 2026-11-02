function refreshPoolWorkerCode(p, repoRoot)
%REFRESHPOOLWORKERCODE Make pool workers resolve GIAnT functions like the client.
%   refreshPoolWorkerCode(p, repoRoot) runs on every worker of the process
%   pool p:
%       addpath(genpath(repoRoot), '-begin'); rehash; clear functions
%   so that (1) this repository shadows any other checkout that defines the
%   same function names, and (2) functions edited since the pool started are
%   re-read from disk.
%
%   Why both are needed
%   -------------------
%   * Path order: a process pool copies the client's path when it STARTS, in
%     the order the client had at that moment, and GIAnT reuses an existing
%     pool whenever its size matches. processSLAP2 puts this repo first on
%     the client path, but a pool left over from another session (for
%     example one started while SLAP2-Voltage-Extraction, which carries its
%     own copies of several GIAnT files, was first) keeps the other order.
%     Workers then silently run the other repo's version of e.g.
%     saveStructToH5 or MultiRoiRegistration>alignAsync.
%   * Function cache: workers keep the compiled copy of a function they have
%     already executed and do not re-check file timestamps the way the
%     client does, so a pool that predates a `git pull` runs the old code.
%
%   Cost is milliseconds per worker, far cheaper than restarting the pool.
%   Thread pools share the client's path and function cache and are skipped.

if nargin < 1 || isempty(p) || ~isvalid(p) || ~isa(p, 'parallel.ProcessPool')
    return
end
if nargin < 2 || isempty(repoRoot)
    repoRoot = fileparts(fileparts(mfilename('fullpath')));   % dependencies/..
end
try
    f = parfevalOnAll(p, @refreshOnWorker, 0, char(repoRoot));
    wait(f);
    if ~isempty(f.Error)
        warning('refreshPoolWorkerCode:Failed', ...
            'Could not refresh worker path/function cache: %s', f.Error.message);
    end
catch ME
    warning('refreshPoolWorkerCode:Failed', ...
        'Could not refresh worker path/function cache: %s', ME.message);
end
end


function refreshOnWorker(repoRoot)
if isfolder(repoRoot)
    addpath(genpath(repoRoot), '-begin');
    rehash;
end
clear functions %#ok<CLFUNC>
end
