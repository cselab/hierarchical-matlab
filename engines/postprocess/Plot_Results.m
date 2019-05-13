%% User parameters

% Assuming dimension > 1 for theta
ind_dim = [1 4; 2 5; 3 6]; % pairs of plotted dimensions [#figure x 2]
f_seed = 1; % index of first figure
size_t = 20; % title font size
size_xy = 16; % x,y label font size
size_line = 2; % line width
size_mark = 10; % marker size
size_dot = 50; % scatter dot base size

% file_name = 'data/tmcmc_single_level';
% file_name = 'data/tmcmc_009';
file_name = 'data/tmcmc_norm_HB';


t_pause = 1; % sec of pause when plotting intermediate samples

%% Plotting final results
load([ file_name '.mat'])

figure(f_seed)
tmp = out_master.N_accept(2:end)./(out_master.N_accept(2:end) + out_master.N_reject(2:end));
plot(2:out_master.gen,tmp,'-xk','LineWidth',size_line,'MarkerSize',size_mark)
title('BASIS statistics','FontSize',size_t)
xlabel('Stage','FontSize',size_xy)
ylabel('Acceptance rate','FontSize',size_xy)
xlim([1 out_master.gen])
ylim([0 1])

for i = 1:size(ind_dim,1)
    figure(f_seed+i)
    clf
    tmp = out_master.lik + out_master.pri;
    scatter(out_master.theta(:,ind_dim(i,1)),out_master.theta(:,ind_dim(i,2)),...
        size_dot*(1+0.1*out_master.Ns),tmp)
    title('Final samples','FontSize',size_t)
    xlabel(['Dim. ',num2str(ind_dim(i,1))],'FontSize',size_xy)
    ylabel(['Dim. ',num2str(ind_dim(i,2))],'FontSize',size_xy)
    
    c = colorbar;
    c1=get(gca,'position');
    c2=get(c,'Position');
    c2(3)=0.5*c2(3);
    set(c,'Position',c2)
    set(gca,'position',c1)
end


%% Plotting intermediate samples

% for k = 1:1
for k = 1:out_master.gen-1
    load([file_name '_gen_'  sprintf('%02d',k)  '.mat'])
    for i = 1:size(ind_dim,1)
        figure(f_seed+i)
        clf
        tmp = runinfo.lik + runinfo.pri;
        scatter(runinfo.theta(:,ind_dim(i,1)),runinfo.theta(:,ind_dim(i,2)),...
            size_dot*(1+0.1*runinfo.Ns),tmp)
        title(['Stage ',num2str(k)],'FontSize',size_t)
        xlabel(['Dim. ',num2str(ind_dim(i,1))],'FontSize',size_xy)
        ylabel(['Dim. ',num2str(ind_dim(i,2))],'FontSize',size_xy)
        
        c = colorbar;
        c1=get(gca,'position');
        c2=get(c,'Position');
        c2(3)=0.5*c2(3);
        set(c,'Position',c2)
        set(gca,'position',c1)
    end
%     pause(t_pause)
    pause
end
%%
load([ file_name '.mat'])

for i = 1:size(ind_dim,1)
    figure(f_seed+i)
    clf
    tmp = out_master.lik + out_master.pri;
    scatter(out_master.theta(:,ind_dim(i,1)),out_master.theta(:,ind_dim(i,2)),...
        size_dot*(1+0.1*out_master.Ns),tmp)
    title('Final samples','FontSize',size_t)
    xlabel(['Dim. ',num2str(ind_dim(i,1))],'FontSize',size_xy)
    ylabel(['Dim. ',num2str(ind_dim(i,2))],'FontSize',size_xy)
    
    c = colorbar;
    c1=get(gca,'position');
    c2=get(c,'Position');
    c2(3)=0.5*c2(3);
    set(c,'Position',c2)
    set(gca,'position',c1)
end