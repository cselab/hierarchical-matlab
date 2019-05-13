function b = start_parallel( N_cpu, spmd_par, spmd_args )

if ~exist('spmd_args','var') || isempty(spmd_args)
  spmd_args = '-W 10:49';
end


b = false;
if( ~spmd_par ) % the machine is not Brutus
    
    if N_cpu > 1
        max_trial = 100;
        tmp = 1;
        tmp_par = gcp('nocreate');
        while (isempty(tmp_par) || (tmp_par.NumWorkers < N_cpu)) && (tmp < max_trial)
            try
                tmp_par = parpool('local', N_cpu);
                pctRunOnAll warning('error', 'MATLAB:ode15s:IntegrationTolNotMet')
                pctRunOnAll warning('off','MATLAB:nearlySingularMatrix')
                b = true;
            catch
                tmp = tmp + 1;
            end
        end
        if tmp == max_trial
            fprintf('matlabpool open failed %d times\n',max_trial)
        end
    end
    
else % the machine is Brutus
    
    if N_cpu > 1
        cluster = parcluster('BrutusLSF8h');
        jobid = getenv('LSB_JOBID');            % If multiple jobs will be run at the same time
        mkdir(jobid);                           % then each one should have its own controlling
        cluster.JobStorageLocation = jobid;     % subdirectory.
        
        %cluster.SubmitArguments = '-W 36:00 -R "rusage[mem=4096]"'; % Workers need less memory.
        cluster.SubmitArguments = spmd_args;
        cluster.NumWorkers = 96;
        
        max_trial = 100;
        tmp = 1;
        tmp_par = gcp('nocreate');
        while (isempty(tmp_par) || (tmp_par.NumWorkers < N_cpu)) && (tmp < max_trial)
            try
                tmp_par = parpool(cluster,N_cpu);
                pctRunOnAll warning('error', 'MATLAB:ode15s:IntegrationTolNotMet' )
                pctRunOnAll warning('off','MATLAB:nearlySingularMatrix')
                b = true;
            catch
                tmp = tmp + 1;
            end
        end
        if tmp == max_trial
            fprintf('parpool open failed %d times on brutus cluster\n',max_trial)
        end
    end
    %==========================================================================
    %==========================================================================
end


end

