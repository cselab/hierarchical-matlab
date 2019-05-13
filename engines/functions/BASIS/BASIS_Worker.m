%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function: local workers at each cpu that calculate likelihood and prior 
%           values in log-scale, perform MCMC chains
% created by Stephen Wu, 2015 April 28
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Note 1: all rejected samples' extra likelihood output not saved
% Note 2: N_accept and N_reject includes burn-in periods

% Last update: 2015 June 9, by Stephen Wu
%   2015-06-02:
%       1. add handling of burn-in period for each MCMC chain [line 73]
%   2015-06-09:
%       1. check work.N_job = 0 [line 25]
%       2. add flexibility on controling burn-in stages [line 71]
%   2015-10-24:
%       1. support len(N_burn)!=Gen_burn [line 73]

function [out_worker] = BASIS_Worker(ind_worker,work,local_para,local_info)

% Initialize empty output
out_worker = struct([]);
if isempty(work) || (work.N_job == 0)
    return
else
    out_worker(work.N_job).theta = [];
end

% Message from worker
if local_para.TF_msg
    fprintf('Worker #%d: My assigned work is %d chains, tot len = %d\n',...
        ind_worker,work.N_job,sum(work.len));
end

if local_info.gen == 1
    % Calculate prior for all sample at once
    tmp = PDF_Pri(work.theta,local_para.pri);
    % Calculate likelihood and store outputs
    for i = 1:work.N_job
        out_worker(i).theta = work.theta(i,:);
        out_worker(i).Ns = 1;
        % Note: out_lik store as a cell here for consistency
        [out_worker(i).lik,out_worker(i).out_lik{1}] = ...
            feval(local_para.lik.name,work.theta(i,:),local_para.lik.para);
        out_worker(i).pri = tmp(i);
        
        % Initialize statistics for acceptance rate
        out_worker(i).N_accept = 0;
        out_worker(i).N_reject = 0;
    end
else
    % Take known likelihood values and do MCMC for each job assigned
    for i = 1:work.N_job
        % initialize the MCMC chain
        thetao = work.theta(i,:);
        lnfo_lik = work.lik(i);
        lnfo_pri = work.pri(i);
        % initialize output from worker (delete later if Ns = 0 at the end)
        out_worker(i).theta = thetao;
        out_worker(i).Ns = 0;
        out_worker(i).lik = lnfo_lik;
        out_worker(i).pri = lnfo_pri;
        out_worker(i).out_lik = work.out_lik(i); % store the cell
        % initialize statistics for acceptance rate
        out_worker(i).N_accept = 0;
        out_worker(i).N_reject = 0;
        
        % Burn-in period
        if local_info.gen <= local_para.Gen_burn
        
        for j = 1:local_para.N_burn(...
                min(local_info.gen,length(local_para.N_burn)))
            % Propose a new sample based on current theta
            [TF_inBDs,thetac] = ...
                Rnd_Prop(thetao,local_info.cov_s,local_para.hard_bds);
            if TF_inBDs
                % Calc. likelihood and prior for proposed sample
                [lnfc_lik,outc_lik] = ...
                    feval(local_para.lik.name,thetac,local_para.lik.para);
                lnfc_pri = PDF_Pri(thetac,local_para.pri);

                % Calc. ratio of acceptance
                if local_para.TF_anneal_pri
                    r = exp(local_info.p(local_info.gen)*...
                        (lnfc_lik + lnfc_pri - lnfo_lik - lnfo_pri));
                else
                    r = exp(local_info.p(local_info.gen)*...
                        (lnfc_lik - lnfo_lik) + (lnfc_pri - lnfo_pri));
                end
                r = min(1,r);
                state = find(mnrnd(1,[r,1-r]));

                % Accept or reject proposed sample
                if state == 1
                    % update thetao
                    thetao = thetac;
                    lnfo_lik = lnfc_lik;
                    lnfo_pri = lnfc_pri;
                    
                    % update initial sample for MCMC seed
                    out_worker(i).theta = thetao;
                    out_worker(i).Ns = 0;
                    out_worker(i).lik = lnfo_lik;
                    out_worker(i).pri = lnfo_pri;
                    out_worker(i).out_lik{1} = outc_lik;
                    
                    out_worker(i).N_accept = out_worker(i).N_accept + 1;
                else
                    out_worker(i).N_reject = out_worker(i).N_reject + 1;
                end
            else
                out_worker(i).N_reject = out_worker(i).N_reject + 1;
            end
        end
        
        end
        
        % Actual MCMC period
        for j = 1:work.len(i)
            % Propose a new sample based on current theta
            [TF_inBDs,thetac] = ...
                Rnd_Prop(thetao,local_info.cov_s,local_para.hard_bds);
            if TF_inBDs
                % Calc. likelihood and prior for proposed sample
                [lnfc_lik,outc_lik] = ...
                    feval(local_para.lik.name,thetac,local_para.lik.para);
                lnfc_pri = PDF_Pri(thetac,local_para.pri);

                % Calc. ratio of acceptance
                if local_para.TF_anneal_pri
                    r = exp(local_info.p(local_info.gen)*...
                        (lnfc_lik + lnfc_pri - lnfo_lik - lnfo_pri));
                else
                    r = exp(local_info.p(local_info.gen)*...
                        (lnfc_lik - lnfo_lik) + (lnfc_pri - lnfo_pri));
                end
                r = min(1,r);
                state = find(mnrnd(1,[r,1-r]));

                % Accept or reject proposed sample
                if state == 1
                    % add new theta to the record & outc_lik
                    out_worker(i).theta(end+1,:) = thetac;
                    out_worker(i).Ns(end+1) = 1;
                    out_worker(i).lik(end+1) = lnfc_lik;
                    out_worker(i).pri(end+1) = lnfc_pri;
                    out_worker(i).out_lik{end+1} = outc_lik;
                    % update thetao
                    thetao = thetac;
                    lnfo_lik = lnfc_lik;
                    lnfo_pri = lnfc_pri;
                    
                    out_worker(i).N_accept = out_worker(i).N_accept + 1;
                else
                    % increase current theta count by one when rejected
                    out_worker(i).Ns(end) = out_worker(i).Ns(end) + 1;
                    
                    out_worker(i).N_reject = out_worker(i).N_reject + 1;
                end
            else
                % increase current theta count by one when out of bounds
                out_worker(i).Ns(end) = out_worker(i).Ns(end) + 1;
                
                out_worker(i).N_reject = out_worker(i).N_reject + 1;
            end
        end
        
        % delete first sample if the count is 0 (from initialization)
        if out_worker(i).Ns(1) == 0
            out_worker(i).theta(1,:) = [];
            out_worker(i).Ns(1) = [];
            out_worker(i).lik(1) = [];
            out_worker(i).pri(1) = [];
            out_worker(i).out_lik(1) = [];
        end
    end 
end

end


%% Function to calculate prior pdf values in log scale
% Note 1: only pass sys_para.pri into para
% Note 2: theta is matrix with size - #samples x #dimension
% Note 3: ln_f(k) = -inf if prior = 0 for the sample (k)
function [ln_f] = PDF_Pri(theta,para)
    if para.TF_custom
        ln_f = feval(para.name_custom_pdf,theta,para.para_custom_pdf);
    else
        ln_f = zeros(size(theta,1),1);
        for i = 1:size(theta,2)
            switch length(para.para{i})
                case 1
                    ln_f = ln_f + ...
                        log(pdf(para.name{i},theta(:,i),para.para{i}));
                case 2
                    ln_f = ln_f + log(pdf(para.name{i},theta(:,i),...
                        para.para{i}(1),para.para{i}(2)));
                case 3
                    ln_f = ln_f + log(pdf(para.name{i},theta(:,i),...
                        para.para{i}(1),para.para{i}(2),para.para{i}(3)));
            end
        end
    end
end


%% Function to propose next sample based on input theta (NOT DONE)
function [TF_inBDs,thetac] = Rnd_Prop(theta,SIG,bds)
    % propose sample with local Gaussian random work
    thetac = mvnrnd(theta,SIG);
    % check hard bounds (assume theta is row vector)
    % bds(1,:) = row vector of min; bds(2,:) = row vector of max
    if any(thetac < bds(1,:)) || any(thetac > bds(2,:))
        TF_inBDs = false;
    else
        TF_inBDs = true;
    end
end

