load ../../data/opt_unif_theta_post_001.mat

figure(1); clf
plotmatrix_hist(out_master.theta(:,:))



load ../../data/HB_unif_theta_post_001.mat

figure(2); clf
plotmatrix_hist(out_master.theta(:,:))