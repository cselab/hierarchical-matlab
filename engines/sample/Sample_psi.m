% Sample hyperparameters psi

% warning off

clear

addpath(genpath('../functions/'));


% uniform prior on theta
sys_para = User_Input_psi_unif;
% gaussian prior on theta
% sys_para = User_Input_psi_norm;


BASIS = BASIS_Master(sys_para);


fprintf( [ repmat('  %f  ',1,size(BASIS.theta,2)) '\n'], mean(BASIS.theta))

