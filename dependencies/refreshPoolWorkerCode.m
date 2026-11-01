function refreshPoolWorkerCode(p)
%REFRESHPOOLWORKERCODE Drop stale compiled functions on existing pool workers.
%   refreshPoolWorkerCode(p) runs CLEAR FUNCTIONS on every worker of the
%   process pool p so the next call re-parses each function from disk.
%
%   Why: the MATLAB client re-checks file timestamps at each prompt, but
%   process-pool workers keep the compiled copy of a function they have
%   already executed. GIAnT reuses an existing pool whenever its size matches,
%   so a pool started before a `git pull` happily runs the OLD code for the
%   whole next run (observed with saveStructToH5: the stack trace named a
%   function that no longer existed in the checkout). Clearing is cheap
%   (milliseconds per worker) compared with restarting the pool.
%
%   Thread pools share the client's function cache and need no refresh.

if nargin < 1 || isempty(p) || ~isvalid(p) || ~isa(p, 'parallel.ProcessPool')
    return
end
try
    f = parfevalOnAll(p, @clearCompiledFunctions, 0);
    wait(f);
    if ~isempty(f.Error)
        warning('refreshPoolWorkerCode:Failed', ...
            'Could not refresh worker function cache: %s', f.Error.message);
    end
catch ME
    warning('refreshPoolWorkerCode:Failed', ...
        'Could not refresh worker function cache: %s', ME.message);
end
end


function clearCompiledFunctions()
clear functions %#ok<CLFUNC>
end
