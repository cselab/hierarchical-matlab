clear; 

data_folder = '../../data/';

k_set = 1;

prefix1 = 'IND_theta_';         prefix2 = 'propagation_IND_';
% prefix1 = 'HB_unif_theta_post_'; prefix2 = 'propagation_post_';

load([ data_folder prefix1 sprintf('%03d',k_set) '.mat']);
load([ data_folder prefix2 sprintf('%03d',k_set) '.mat']);

th = out_master.theta; 
clear out_master;



n1 = size(y,1);
n2 = size(y,2);


sigma = th(1:n1,end);



%%

v = [ 0.05, 0.95, 0.005, 0.995] ;

p = zeros(length(v),n2);


for i = 1 : n2

    mu = y(:,i);
    f = @(x) sum( normcdf(x,mu,sigma) )/n1;    
    
    
    for k=1:length(v)
        p(k,i) = fzero(@(x)f(x)-v(k),0.5);
    end
    
end

%% PLOT confidence intervals


load([ data_folder 'data_set_' sprintf('%03d',k_set) '.mat']);


fg = figure(1); clf

fg.Position = [417   667   853   678];

cols = [    0         0.4470    0.7410
            0.8500    0.3250    0.0980
            0.9290    0.6940    0.1250
            0.4940    0.1840    0.5560
            0.4660    0.6740    0.1880
            0.3010    0.7450    0.9330
            0.6350    0.0780    0.1840];

mf = mean(y);



% 
bounds(:,:,2) = [p(1,:)-mf; mf-p(2,:)]';
bounds(:,:,1) = [p(3,:)-mf; mf-p(4,:)]';
[ll,pp] = boundedline( x, [mf;mf], bounds, 'cmap', cols); hold on


pd = plot( data.x, data.y, 'k*' ); 
grid on;


lgnd = legend( [pd ll(2) pp(2) pp(1)], '  measurements', '  expected prediction' , '  95% credible interval', '  99% credible interval' );
lgnd.Location = 'southeast';


grid on

ax = gca;
axis tight