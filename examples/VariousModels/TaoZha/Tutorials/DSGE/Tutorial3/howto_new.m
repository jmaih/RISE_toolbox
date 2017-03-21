%%
% RISE Tutorial by Dr. Tao Zha
%% housekeeping
clear
close all
clc
% This tutorial needs:
% - report
%% Instructions
% - Please run this file block by block and make sure you read the
% comments in each block to understand what it does. If there is anything
% you do not understand, ask questions to the instructor or to your
% neighbor
% - to run a particular block, click on the block and then on your keyboard
% press CTRL+Enter
%% add the paths to RISE, the data and the models
setpaths=true;
if setpaths
    addpath([pwd,filesep,'Models']) % folder with the models
    addpath([pwd,filesep,'Data']) % folder containing the data
    % Instead of copying the folder in github somewhere on your disk, you
    % can alternatively just set the paths so that matlab reads information
    % from the examples folder but does not write to them. Here is an
    % example of how to do that. It assumes that your main folder is not
    % under github!
end
%% Bring in some data and transform them into rise_time_series

tmp=load('Data/data_nk3eq_8501_1301');  %qdatae
dataList={
    'Y','output gap  y_t (log GDP_t - log GDPPotential_t)'
    'PAI','PCE core inflation pi_t (log P_t - log P_{t-1})'
    'I','FFR R_t log(1+ffr/400)  (quarterly rate already)'
    };
mydata=struct();
startdate='1985Q1';
for id=1:size(dataList,1)
    % we just give the start date, RISE automatically understands that we
    % are dealing with quarterly data by the format startdate
    mydata.(dataList{id,1})=ts(startdate,... start date
        tmp.qdatae(:,id+1),... the data
        dataList{id,2}); % list of variables
end

%% plot your data, compute basic statistics and look at both carefully
varlist=fieldnames(mydata);
figure('name','US data')
nvars=numel(varlist);
for id=1:nvars
    subplot(nvars,1,id)
    dd=mydata.(varlist{id});
    plot(dd,'linewidth',2)
    title(mydata.(varlist{id}).varnames)
    fprintf('%s:: mean %0.3f  stdev %0.3f\n',mydata.(varlist{id}).varnames{1},mean(dd),std(dd));
end
[~,tmp]=sup_label(['US data ',mydata.(varlist{id}).start,':', mydata.(varlist{id}).finish],'t');
set(tmp,'fontsize',15)
%% Read the model(s)

model_names={'svar_constant','svar_policy','svar_volatility',...
    'svar_policy_volatility'};
nmodels=numel(model_names);
estim_models=cell(1,nmodels);

% rather than putting all the models in the same vector as we did earlier,
% we put them in a cell array. If we put them in the same vector and call
% the estimation function, RISE will think that we want to estimate a
% pareto-type of model. But this is not what we want to do and is probably
% beyond the scope of these lectures.

% we loop through the different models using the information in the labels
for imod=1:nmodels
    % replace "for" by "parfor" if you want to use parallel computation
    estim_models{imod}=rise(model_names{imod},... % name of the file to read
        'saveas',true,... % we ask rise to write the expanded model to disk
        'data',mydata,... % we may assign the data now or later
        'solve_linear',true... take advantage of the fact that the models are conditionally linear
        );
    % a model with multiple files inserted can be difficult to read. The
    % expanded model could be useful for understanding what RISE does and
    % for debugging purposes. The expanded model contains all the details
    % of the individual files (without the comments)
end

%% We estimate the models or filter them directly
close all
% if we have the parallel computing toolbox, we can estimate all models in
% one go
% parpool(nmodels)
% parfor imod=1:nmodels %
parfor imod=1:nmodels %
    % replace "for" by "parfor" if you want to use parallel computation
    disp('*--------------------------------------------------------------*')
    disp(['*-----Estimation of ',model_names{imod},' model-------*'])
    disp('*--------------------------------------------------------------*')
    estim_models{imod}=estimate(estim_models{imod},'optimizer','fmincon');
    close all
end


%% Simulation of the posterior distribution
% obj below is a rise object with the simulation statistics included in
% obj{imod}.estimation.posterior_simulation for imod=1:nmodels
clc
Results=cell(1,nmodels);
objective=cell(1,nmodels);
lb=cell(1,nmodels);
ub=cell(1,nmodels);
ndraws_mcmc         = 1500;  % number of parameter draws through MCMC.
ndraws_burnin       = floor(0.1*ndraws_mcmc); % number of parameter draws to be burned
mcmc_options=struct('burnin',ndraws_burnin,'N',ndraws_mcmc,'thin',1);

parfor imod=1:nmodels
    % replace "for" by "parfor" if you want to use parallel computation
    disp('*--------------------------------------------------------------*')
    disp(['*----------- MCMC for ',model_names{imod},' model-------------*'])
    disp('*--------------------------------------------------------------*')
    % Here the simulations are stored in variable, which here we call
    % postSims.
    %-----------------------------------------------------------------
    [objective{imod},lb{imod},ub{imod},x0,SIG]=pull_objective(estim_models{imod});
    
    Results{imod}=mh_sampler(objective{imod},lb{imod},ub{imod},mcmc_options,x0,SIG);
    
end
%% Marginal Likelihood
lmdd=cell(1,nmodels);
do_log_marginal_data_density=~false;
if do_log_marginal_data_density
    for imod=1:nmodels
        % replace "for" by "parfor" if you want to use parallel computation
        disp('*--------------------------------------------------------------*')
        disp(['*--- Marginal data density for ',model_names{imod},' model---*'])
        disp('*--------------------------------------------------------------*')
        
        lmdd{imod} = mcmc_mdd(Results{imod}.pop,lb{imod},ub{imod},...
            struct('log_post_kern',objective{imod},... % function to MINIMIZE !!!
            'algorithm','mhm'));% MDD algorithm
        
    end
end


