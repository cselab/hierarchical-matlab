clear;

load ../../data/CMA_psi.mat

psi = bestever.x;

Ns = 10000;

theta = [ unifrnd(psi(1),psi(2),Ns,1) , unifrnd(psi(3),psi(4),Ns,1), unifrnd(psi(5),psi(6),Ns,1), unifrnd(psi(7),psi(8),Ns,1)    ];

save('../../data/theta_opt.mat','theta')