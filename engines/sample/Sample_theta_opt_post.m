% Main file to run BASIS

% warning off


clear

addpath('../functions/BASIS/');

% Nsets = 1;

tic

for k = 1


    sys_para = User_Input_unif_theta_opt_post( k );
    BASIS = BASIS_Master(sys_para);

    fprintf( [ repmat('  %f  ',1,size(BASIS.theta,2)) '\n'], mean(BASIS.theta))
    
end

toc
