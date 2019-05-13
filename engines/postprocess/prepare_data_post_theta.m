clear;

index = 1;
s_index = sprintf('%03d',index);

data_folder = '../../data/';


theta_file_name = [ data_folder 'IND_theta_' s_index '.mat'];
psi_file_name   = [ data_folder 'HB_unif_psi.mat'];

th_samples  = load(theta_file_name);
psi_samples = load(  psi_file_name);

ind_evidence = exp(th_samples.out_master.lnEv);

Ns_theta = size( th_samples.out_master.theta,1);
Ns_psi   = size(psi_samples.out_master.theta,1);

B = exp(th_samples.out_master.pri);

denom = zeros(Ns_psi,1);

for k = 1:Ns_psi
   disp(k)
    psi = psi_samples.out_master.theta(k,:);
    
    
    prior_theta_psi = @(x)  unifpdf( x(:,1), psi(1), psi(1) + psi(2))...
                         .* unifpdf( x(:,2), psi(3), psi(3) + psi(4))...
                         .* unifpdf( x(:,3), psi(5), psi(5) + psi(6))...
                         .* unifpdf( x(:,4), psi(7), psi(7) + psi(8));

    theta = th_samples.out_master.theta;
    A = prior_theta_psi(theta);
    
    
    denom(k) = sum(A./B);
    
end



save([data_folder 'aux_' s_index '.mat'],'denom')
