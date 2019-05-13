%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function: main controller of BASIS algorithm that get initial samples, 
%           pack and distribute work, collect and analyze results
% created by Stephen Wu, 2015 April 28
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Last update: 2015 Oct 24, by Stephen Wu
%   2015-05-04:
%       1. update SPMD code - output(i) -> output(i,1) [line 204]
%   2015-05-19:
%       1. update len_chain calc. to reduce storage needed [line 151,154]
%   2015-05-27:
%       1. update Work_Balance to have uniform chain length [line 221]
%   2015-06-09:
%       1. check work.N_job = 0, add warning [line 176]
%   2015-06-17:
%       1. add run time tracking [line 27, 102]
%   2015-09-18:
%       1. add min. ini. value for optimization of p [line 348]
%   2015-10-24:
%       1. support len(sys_para.maxP/minP)!=N_gen_maxP/minP [line 326,332]

function [out_master] = BASIS_Master(sys_para)

tmp_t = tic;

%% Initialization
runinfo.gen = 0; % initialize as 0, will consider prior samples as stage 1
runinfo.p = 0; % initialize stage 1 power = 0
runinfo.TF_fin = false; % true when p -> 1 (last stage)
runinfo.cov_s = []; % initialize for stage 1

%% Generate initial samples
if sys_para.TF_prop_pri
    runinfo.theta = sys_para.prop_theta;
elseif sys_para.pri.TF_custom
    runinfo.theta = feval(sys_para.pri.name_custom_rnd,...
        sys_para.pri.para_custom_rnd);
else
    runinfo.theta = zeros(sys_para.N_s,sys_para.N_dim);
    for i = 1:sys_para.N_dim
        switch length(sys_para.pri.para{i})
            case 1
                runinfo.theta(:,i) = random(sys_para.pri.name{i},...
                    sys_para.pri.para{i},sys_para.N_s,1);
            case 2
                runinfo.theta(:,i) = random(sys_para.pri.name{i},...
                    sys_para.pri.para{i}(1),sys_para.pri.para{i}(2),...
                    sys_para.N_s,1);
            case 3
                runinfo.theta(:,i) = random(sys_para.pri.name{i},...
                    sys_para.pri.para{i}(1),sys_para.pri.para{i}(2),...
                    sys_para.pri.para{i}(3),sys_para.N_s,1);
        end
    end
end


%% Main loop
while runinfo.gen < sys_para.max_gen
    % increment to current stage
    runinfo.gen = runinfo.gen + 1;
    
    % pack and distribute work to local workers, then get local results
    [out_worker] = Pack_And_Ship_Work(runinfo,sys_para);
    % assemble results from local workers
    [runinfo] = Collect_Work(out_worker,runinfo,sys_para);
    % update annealing power and evidence (not needed for final stage)
    if runinfo.TF_fin
        % quit the loop when p = 1
        runinfo.lnEv = sum(runinfo.S_lnEv);
        break
    else
        [runinfo] = Calc_Stat(runinfo,sys_para);
    end
    
    % messaging and intermediate file saving
    if sys_para.TF_msg
        fprintf('Master: Stage %d completed (p = %.2e)\n',runinfo.gen,...
            runinfo.p(runinfo.gen));
    end
    if sys_para.TF_save
        save([sys_para.save_file '_gen_',sprintf('%02d',runinfo.gen),'.mat'],'runinfo')
        if sys_para.TF_msg
            fprintf('   Stage %d saved\n',runinfo.gen);
        end
    end
end

% save final results
out_master = runinfo;
    % get mean and cov of the final samples
tmp = runinfo.Ns/sum(runinfo.Ns);
out_master.mean_fin = tmp' * runinfo.theta;
out_master.cov_fin = runinfo.theta - ...
    repmat(tmp'*runinfo.theta,size(runinfo.theta,1),1); 
out_master.cov_fin = out_master.cov_fin' * (out_master.cov_fin .* ...
    repmat(tmp,1,sys_para.N_dim));
out_master.cov_fin = 0.5 * (out_master.cov_fin + out_master.cov_fin');

out_master.runtime = toc(tmp_t);

save_file = [sys_para.save_file '.mat'];

save( save_file ,'out_master','sys_para');
if sys_para.TF_msg
    if runinfo.TF_fin
        fprintf('Master: Final stage finished successfully\n');
    else
        fprintf('Master: Max. stages reached, program exited\n');
    end
    fprintf('   Final results saved in %s\n',save_file);
end

end


%% Function to pack and ship work
function [out_worker] = Pack_And_Ship_Work(runinfo,sys_para)

% initialize 'work' as N_cpu x 1 empty struct
work = struct;
work(sys_para.N_cpu).N_job = 0;

% Distribute work
if runinfo.gen == 1
    % equally distribute work for likelihood and prior calculation
    
    % record number of jobs per cpu equally
    N_job = zeros(sys_para.N_cpu,1);
    N_job(1:mod(sys_para.N_s,sys_para.N_cpu)) = 1;
    N_job = N_job + floor(sys_para.N_s/sys_para.N_cpu);
    
    % distribution work based on N_job (only theta needed at gen = 1)
    tmp = 0;
    for i = 1:sys_para.N_cpu
        work(i).N_job = N_job(i);
        work(i).theta = runinfo.theta(tmp+1:tmp+N_job(i),:);
        work(i).len = ones(N_job(i),1);
        
        tmp = tmp + N_job(i);
    end
    
    if sys_para.TF_msg
        fprintf('Master: Stage 1 work packed\n')
    end
    
else
    % resample to determine length of each chain
        % len_chain: column vector
    if runinfo.TF_fin
        len_chain = mnrnd(sys_para.N_s_fin,runinfo.w)';
%         len_chain = sum(mnrnd(1,runinfo.w,sys_para.N_s_fin))';
    else
        len_chain = mnrnd(sys_para.N_s,runinfo.w)';
%         len_chain = sum(mnrnd(1,runinfo.w,sys_para.N_s))';
    end
    
    % distribute work smartly!
    [job_ind,job_len,stat] = ...
        Work_Balancing(len_chain,sys_para.N_cpu,sys_para.max_cLen);
    
    % pack the work
    for i = 1:sys_para.N_cpu
        work(i).N_job = length(job_ind{i});
        work(i).theta = runinfo.theta(job_ind{i},:);
        work(i).len = job_len{i};
        work(i).lik = runinfo.lik(job_ind{i},:);
        work(i).pri = runinfo.pri(job_ind{i},:);
        work(i).out_lik = runinfo.out_lik(job_ind{i});
    end
    
    if sys_para.TF_msg
        fprintf('Master: Stage %d work packed\n',runinfo.gen);
        fprintf('   Mean job length = %f, std = %f\n',stat.mean,stat.std);
        tmp = [];
        for i = 1:sys_para.N_cpu
            if work(i).N_job == 0
                tmp(end+1) = i;
            end
        end
        if ~isempty(tmp)
            fprintf('   Warning: worker');
            fprintf(' %d',tmp);
            fprintf(' has 0 jobs\n');
        end
    end
end

out_worker = cell(sys_para.N_cpu,1);
% copy neccessary runinfo to local_info (reduce data transfer burden)
local_info.gen = runinfo.gen;
local_info.p = runinfo.p;
local_info.cov_s = runinfo.cov_s;
if sys_para.N_cpu > 1
    parfor i = 1:sys_para.N_cpu
        out_worker{i} = BASIS_Worker(i,work(i),sys_para,local_info);
    end
    
%     % SPMD style parallelization
%     spmd
%         out_worker = cell(sys_para.N_cpu,1,codistributor());
%         for i = drange(1:sys_para.N_cpu)
%            % important to be (i,1) not (i) - or else spmd won't work!
%            out_worker(i,1)={BASIS_Worker(i,work(i),sys_para,local_info)};
%         end
%     end %of spmd
%     out_worker = gather(out_worker);
    
else
    for i = 1:sys_para.N_cpu
        out_worker{i} = BASIS_Worker(i,work(i),sys_para,local_info);
    end
end

end


%% Function to distribute work equally among available CPUs
function [job_ind,job_len,stat] = Work_Balancing(len_chain,N_cpu,max_cLen)
    % break chains that exceed max. chain length
    ind_bad = find(len_chain > max_cLen);
    ind_cLen = [1:length(len_chain)]'; % use column vector
    ind_cLen(ind_bad) = [];
    cLen_full = len_chain;
    cLen_full(ind_bad) = [];
    
    % New version: evenly distributed chains
    for i = 1:length(ind_bad)
        tmp = floor(len_chain(ind_bad(i))/max_cLen);
        if mod(len_chain(ind_bad(i)),max_cLen) ~= 0
            tmp = tmp + 1;
        end
        len = length(cLen_full);
        cLen_full(len+1:len+tmp) = floor(len_chain(ind_bad(i))/tmp);
        ind_cLen(len+1:len+tmp) = ind_bad(i);
        % add remaining length evenly
        tmp = mod(len_chain(ind_bad(i)),tmp);
        if tmp ~= 0
            cLen_full(len+1:len+tmp) = cLen_full(len+1:len+tmp) + 1;
        end
    end
    
%     % Old version: non-evenly distributed chains
%     for i = 1:length(ind_bad)
%         tmp = floor(len_chain(ind_bad(i))/max_cLen);
%         cLen_full(end+1:end+tmp) = max_cLen;
%         ind_cLen(end+1:end+tmp) = ind_bad(i);
%         tmp = mod(len_chain(ind_bad(i)),max_cLen);
%         if tmp ~= 0
%             cLen_full(end+1) = tmp;
%             ind_cLen(end+1) = ind_bad(i);
%         end
%     end
    
    % distribute work equally
    job_ind = cell(N_cpu,1);
    job_len = cell(N_cpu,1);
    cpu_load = zeros(N_cpu,1);
    [cLen_sort,tmp] = sort(cLen_full,'descend'); % zeros at the end
    ind_cLen = ind_cLen(tmp);
        % distribute only non-zero length chains
    for i = 1:sum(cLen_sort > 0)
        [~,tmp] = min(cpu_load);
        cpu_load(tmp) = cpu_load(tmp) + cLen_sort(i);
        job_ind{tmp}(end+1) = ind_cLen(i);
        job_len{tmp}(end+1) = cLen_sort(i);
    end
    
    % record statistics of the work distribution
    stat.mean = mean(cpu_load);
    stat.std = std(cpu_load);
end


%% Function to assemble local results
function [runinfo] = Collect_Work(out_worker,runinfo,sys_para)

% count number of total distinct sample
tmp = 0;
for i = 1:sys_para.N_cpu
    for j = 1:length(out_worker{i})
        tmp = tmp + length(out_worker{i}(j).Ns);
    end
end

% re-initialize part of runinfo
runinfo.theta = zeros(tmp,sys_para.N_dim);
runinfo.Ns = zeros(tmp,1);
runinfo.pri = zeros(tmp,1);
runinfo.lik = zeros(tmp,1);
runinfo.out_lik = cell(tmp,1);
runinfo.N_accept(runinfo.gen) = 0;
runinfo.N_reject(runinfo.gen) = 0;

% assemble output from workers
tmp = 0;
for i = 1:sys_para.N_cpu
    for j = 1:length(out_worker{i})
        len_out = length(out_worker{i}(j).Ns);
        runinfo.theta(tmp+1:tmp+len_out,:) = out_worker{i}(j).theta;
        runinfo.Ns(tmp+1:tmp+len_out) = out_worker{i}(j).Ns;
        runinfo.pri(tmp+1:tmp+len_out) = out_worker{i}(j).pri;
        runinfo.lik(tmp+1:tmp+len_out) = out_worker{i}(j).lik;
        runinfo.out_lik(tmp+1:tmp+len_out) = out_worker{i}(j).out_lik;
        runinfo.N_accept(runinfo.gen) = ...
            runinfo.N_accept(runinfo.gen) + out_worker{i}(j).N_accept;
        runinfo.N_reject(runinfo.gen) = ...
            runinfo.N_reject(runinfo.gen) + out_worker{i}(j).N_reject;

        tmp = tmp + len_out;
    end
end

end


%% Function to calculate statistics and ready for next stage
function [runinfo] = Calc_Stat(runinfo,sys_para)

% control lower and upper limit of the annealing power with fixed
% number of stages if specified, i.e., value ~= 0 
% (Note: for prior stage, runinfo.gen = 1)
if logical(sys_para.N_gen_minP) && ...
        (runinfo.gen <= sys_para.N_gen_minP)
    tmp_min = sys_para.minP(min(runinfo.gen,length(sys_para.minP)));
else
    tmp_min = 0;
end
if logical(sys_para.N_gen_maxP) && ...
        (runinfo.gen <= sys_para.N_gen_maxP)
    tmp_max = sys_para.maxP(min(runinfo.gen,length(sys_para.maxP)));
else
    tmp_max = 1;
end

% optimization on next power based on c.o.v. of samples
if sys_para.TF_anneal_pri
    lnf = runinfo.lik + runinfo.pri;
else
    lnf = runinfo.lik;
end
    % include repeated samples
tmp = find(runinfo.Ns > 1);
for i = 1:length(tmp)
    lnf(end+1:end+runinfo.Ns(tmp(i))-1) = lnf(tmp(i));
end
   % optimization with fmincon 
[runinfo.p(runinfo.gen+1),runinfo.cv_err2(runinfo.gen)] = fmincon(@(x) ...
    Objlogp(x,lnf,runinfo.p(runinfo.gen),sys_para.cv_tol),...
    runinfo.p(runinfo.gen)+tmp_min+sys_para.opt_iniP,[],[],[],[],...
    runinfo.p(runinfo.gen)+tmp_min,runinfo.p(runinfo.gen)+tmp_max,...
    [],sys_para.opt_setup);

% check if next stage is the last stage and set p = 1
if runinfo.p(runinfo.gen+1) >= 1
    runinfo.TF_fin = true;
    runinfo.p(runinfo.gen+1) = 1;
end

% update weights and calculate partial ln(evidence)
    % case dependent formula for the ratio of weights
if sys_para.TF_anneal_pri
    if runinfo.gen == 1
        if sys_para.TF_prop_pri
            tmp = max( (runinfo.lik+runinfo.pri)*...
                runinfo.p(runinfo.gen+1) - sys_para.prop_val);
            runinfo.w = exp( (runinfo.lik+runinfo.pri)*...
                runinfo.p(runinfo.gen+1) - sys_para.prop_val - tmp );
        else
            tmp = max( (runinfo.lik+runinfo.pri)*...
                runinfo.p(runinfo.gen+1) - runinfo.pri);
            runinfo.w = exp( (runinfo.lik+runinfo.pri)*...
                runinfo.p(runinfo.gen+1) - runinfo.pri - tmp );
        end
    else
        tmp = max( (runinfo.lik+runinfo.pri)*(runinfo.p(runinfo.gen+1)...
            - runinfo.p(runinfo.gen)) );
        runinfo.w = exp( (runinfo.lik+runinfo.pri)*...
            (runinfo.p(runinfo.gen+1)-runinfo.p(runinfo.gen)) - tmp );
    end
else
    if (runinfo.gen == 1) && sys_para.TF_prop_pri
        tmp = max( runinfo.lik*runinfo.p(runinfo.gen+1) + runinfo.pri...
            - sys_para.prop_val);
        runinfo.w = exp( runinfo.lik*runinfo.p(runinfo.gen+1) + ...
            runinfo.pri - sys_para.prop_val - tmp );
    else
        tmp = max( runinfo.lik*(runinfo.p(runinfo.gen+1)...
            - runinfo.p(runinfo.gen)) );
        runinfo.w = exp( runinfo.lik*...
            (runinfo.p(runinfo.gen+1)-runinfo.p(runinfo.gen)) - tmp );
    end
end
% scale weights with number of repeated samples
runinfo.w = runinfo.w .* runinfo.Ns;
% calculate partial ln(evidence)
runinfo.S_lnEv(runinfo.gen) = log(sum(runinfo.w))+tmp-log(sys_para.N_s);
% normalize w
runinfo.w = runinfo.w / sum(runinfo.w);

% calculate covariance matrix for MCMC proposals in next stage
    % remove weighted mean
runinfo.cov_s = runinfo.theta - ...
    repmat(runinfo.w'*runinfo.theta,size(runinfo.theta,1),1); 
    % weighted sample cov
runinfo.cov_s = runinfo.cov_s' * (runinfo.cov_s .* ...
    repmat(runinfo.w,1,sys_para.N_dim)); 
    % guarantee symmetry of cov. matrix & scale with beta^2
runinfo.cov_s = (sys_para.beta2 * 0.5) * ...
    (runinfo.cov_s + runinfo.cov_s'); 
    
end


%% Function to calculate cov given sample likelihoods
function [CoefVar] = Objlogp(x,fj,pj,tol)
    fjmax = max(fj);
    q = exp( (fj-fjmax)*(x-pj) );
    q = q / sum(q);
    % CoefVar = abs(std(q)/mean(q) - tol);
    CoefVar = (std(q)/mean(q) - tol)^2;  
end

