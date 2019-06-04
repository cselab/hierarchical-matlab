% sample the parameters for all the data set assuming independency, 
% i.e., sample from p(theta_i|data_i,model_i)

warning('off','MATLAB:dispatcher:UnresolvedFunctionHandle')

clear

addpath('../functions/BASIS/');

Nsets =   10;

tic

for k=1:Nsets

    fprintf('Bayes on data set: %03d / %4d \n',k,Nsets)

    sys_para = User_Input_theta_ind( k );
    BASIS = BASIS_Master(sys_para);

    fprintf( [ repmat('  %f  ',1,size(BASIS.theta,2)) '\n'], mean(BASIS.theta))
    
end

toc
