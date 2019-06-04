% infer the parameters for the k-th data set by optimizing the likelihood
clc; clear

warning('off','MATLAB:dispatcher:UnresolvedFunctionHandle')


k = 1;


data_folder = '../../data/';


addpath('../functions/tools')
addpath('../functions/CMA')
addpath('../functions/logistic_model/')

loglike_func = 'loglike';


file_name = [ data_folder 'data_set_' sprintf('%03d',k) '.mat'];

data = load(file_name);

data.CMA = true;



%% Options
% change warning of ode15s on Tolerance to error for catching
warnId = 'MATLAB:ode15s:IntegrationTolNotMet';
warnstate = warning('error', warnId);

start_parallel(2, false );

opts.CMA.active = 0;
opts.PopSize = 100; % population size in every generation
opts.Resume = 0;
opts.MaxFunEvals = 300000;
opts.LBounds = [ 200,    0,  0,  1e-3 ]';    % upper bound in search space
opts.UBounds = [ 400,  100,  5,  50]';       % lower bound in search space
opts.Noise.on = 0;
opts.LogModulo = 1;
opts.LogPlot = 1;
opts.DispModulo = 1;
opts.EvalParallel = 1;
opts.EvalInitialX = 1;
opts.TolX = 1e-8;


opts.SaveFilename = [ data_folder 'opt_theta' sprintf('%03d',k) '.mat'];
opts.LogFilenamePrefix = [ data_folder 'outcmaes_'];


xinit = [ 300 , 50, 0, 1]';


CMA = cmaes_parfor( loglike_func,  xinit,[], opts,data);
%%
disp(CMA)
% end

