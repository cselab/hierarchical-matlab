# hierarchical-matlab
Matlab implementation of Hierarchical Bayesian Inference


Steps for a complete Hierarchical Bayesian sampling with propagation of
uncertainty in the predictions


step 0:
    make synthetic data.
    go to engines/ and run make_data.m



step 0.1: (optional)
    optimize the parameters for each of the data sets by running CMA
    go to /optimize/
    run optimize_theta_ind.m
    compare the inferred variables in CMA with the nominal in data.data.theta and data.data.std_data


step 1:
    sample individual parameters using the TMCMC algorithm
    go to /sample/ and run Sample_theta_ind.m  (it will take some time)
    go to /postprocess/ and run "load ../../data/IND_theta_002.mat; plotmatrix_hist(out_master.theta);" to see a histogram of the parameters for the 2nd data set


step 2:
    sample the hyperparameters psi using uniform prior for psi
	go to /sample/ and run Sample_psi.m
    go to /postprocess/ and run "load ../../data/HB_unif_psi.mat; plotmatrix_hist(out_master.theta);" to see a histogram of psi


step 3a:
    run the postprocessing tool
    go to /postprocess/ and run prepare_data_post_theta.m


step 3b:
    sample the posterior distribution of theta_i
    go to /sample/ and run Sample_theta_post.m

step 4:
    estimate the prediction
    go to /propagate/
    run run_model.m to run the model with theta samples
    run uq_bounds to plot the uncertainty in the prediction
