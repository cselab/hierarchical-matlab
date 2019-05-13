clear
addpath('../functions/logistic_model/');
addpath(genpath('../functions/plot_tools/'));

index = 3;
s_index = sprintf('%03d',index);

data_folder = '../data/';
load([data_folder 'HB_unif_theta_post_' s_index '.mat']);
load([data_folder 'data_set_' s_index '.mat']);

save_file = [data_folder 'HB_quantiles_' s_index '.mat'];
LOAD = false;

%%
Ns = 1000;%size(out_master.theta,1);

t = data.x;

f = zeros(Ns,length(t));

for i=1:Ns
    f(i,:) = my_model(t,out_master.theta(i,:),1);
end



sigma = out_master.theta(1:Ns,end)';
tmp = zeros(1,1,Ns);
tmp(1,1,:) = sigma;
sigma = tmp;
p = ones(1,Ns)/Ns;

for i=1:length(t)
    obj{i} = gmdistribution( f(:,i),sigma,p);
end






%% ========================================================================
i=1;


v1 = 0.05; v2 = 0.95; v3 = 0.01; v4 = 0.99;

if( LOAD )
    load(save_file);
else
    for i=1:length(t)

        [i length(t)]

        p1(i) = fzero(@(x)obj{i}.cdf(x)-v1,0.5);
        p2(i) = fzero(@(x)obj{i}.cdf(x)-v2,0.5);

        p3(i) = fzero(@(x)obj{i}.cdf(x)-v3,0.5);
        p4(i) = fzero(@(x)obj{i}.cdf(x)-v4,0.5);
        
    end
    save(save_file,'p1','p2','p3','p4');
end

%%
figure(1); clf

x=10:0.01:60;
plot(x,obj{1}.cdf(x'));  hold on

plot([p1(1) p1(1)],[0 obj{1}.cdf(p1(1))],'r')
plot([x(1) p1(1)],[v1 v1],'r')

plot([p2(1) p2(1)],[0 obj{1}.cdf(p2(1))],'r');
plot([x(1) p2(1)],[v2 v2],'r');




axis tight
%% ========================================================================
figure(2); clf
plot(x,obj{1}.pdf(x')); hold on
plot([p1(1) p1(1)],[0 obj{1}.pdf(p1(1))],'r'); 
plot([p2(1) p2(1)],[0 obj{1}.pdf(p2(1))],'r')

axis tight
%% ========================================================================

%% ========================================================================
figure(4); clf

cols = [         0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];

mf = mean(f);

bounds(:,:,1) = [p1-mf; mf-p2]';
bounds(:,:,2) = [p3-mf; mf-p4]';
[l,p] = boundedline(t, [mf;mf], bounds, 'cmap', cols([7 3],:) , 'alpha'); hold on

outlinebounds(l,p);

plot(data.x,data.y,'ko'); 
grid on