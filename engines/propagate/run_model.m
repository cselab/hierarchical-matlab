
clear

data_folder = '../../data/';

k = 1;

prefix1 = 'IND_theta_';         prefix2 = 'propagation_IND_';
% prefix1 = 'HB_unif_theta_post_'; prefix2 = 'propagation_post_';

load([ data_folder prefix1 sprintf('%03d',k) '.mat']);
load([ data_folder 'data_set_' sprintf('%03d',k) '.mat']);

Ns = 5000;

T = data.x(end);

x = 0:0.01:T+2;

y = zeros(Ns,length(x));

for i=1:Ns
    
    y(i,:) = data.modelfun(x, out_master.theta(i,1:data.Np));
    
end
    


save_file = [ data_folder prefix2 sprintf('%03d',k) '.mat'];

save( save_file, 'x', 'y'  )