% Main file to run BASIS

% warning off


clear

addpath('../functions/BASIS/');

% Nsets = 1;

tic

for k = 1:5

%     fprintf('Bayes on data set: %03d / %4d \n',k,Nsets)

    sys_para = User_Input_unif_theta_post( k );
    BASIS = BASIS_Master(sys_para);

    fprintf( [ repmat('  %f  ',1,size(BASIS.theta,2)) '\n'], mean(BASIS.theta))
    
end

toc
% delete(gcp)
% exit