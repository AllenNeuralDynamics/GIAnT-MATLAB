function pool = extractionPool(params)
%EXTRACTIONPOOL a parallel pool of worker processes for source extraction,
%reusing a process pool that is already running (e.g. SILo's) if it is
%large enough
%
%   Subproblems scale better on processes than on threads, which contend
%   inside one MATLAB process: with 16 subproblems running at once, each ran
%   2.6-3.8x slower than alone on processes, and 5-6x slower on threads.
%   Starting 16 process workers on a 16-core Windows workstation often
%   failed (a worker exits during startup), while 12 started reliably and
%   extracted as fast. The pool has params.nWorkers workers (12 by default),
%   at most the 'Processes' profile's. A running process pool with at least
%   that many is reused; a smaller one (SILo limits its pool by the memory
%   its trials need, which subproblems do not) is replaced. A failed start
%   is retried; if none succeeds, a thread pool is used.
nWorkers = 12;
if nargin && isfield(params, 'nWorkers') && ~isempty(params.nWorkers)
    nWorkers = params.nWorkers;
end
nWorkers = max(1, min(nWorkers, parcluster('Processes').NumWorkers));
pool = gcp('nocreate');
if isa(pool, 'parallel.ProcessPool') && pool.NumWorkers >= nWorkers
    return
end
delete(pool);
for attempt = 1:3
    try
        if attempt > 1
            pause(30); %a failed start can report its workers' exit late; let it, inside this try
        end
        pool = parpool('Processes', nWorkers);
        return
    catch err
        delete(gcp('nocreate'));
        fprintf(2, 'extractionPool: starting %d worker processes failed (%s)\n', nWorkers, err.message);
    end
end
warning('extractionPool:threads', 'Could not start a process pool; using threads');
pool = parpool('Threads');
end
