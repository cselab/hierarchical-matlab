% plots the data created with the make_data script
%
%

clear
data_folder = '../../data/';

warning('off', 'MATLAB:dispatcher:UnresolvedFunctionHandle');

figure(1); clf

N = 5; % number of data files

for k=1:N
    
    load([ data_folder 'data_set_' sprintf('%03d',k) '.mat']);
    
    plot(data.x,data.y,'-o'); hold on
    
end
grid on

ax=gca;

ax.XLabel.String = 'time';
ax.YLabel.String = 'Population Size';

ax.FontSize = 20;