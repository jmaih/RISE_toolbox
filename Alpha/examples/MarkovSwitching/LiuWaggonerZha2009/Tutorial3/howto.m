%%
% RISE Tutorial by Dr. Tao Zha
%% housekeeping
clear all
close all
clc
% This tutorial needs:
% - updating 
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
    % change th path below according to your own system
%     addpath('C:\Users\L5\Documents\GitHub\RISE_toolbox')
%     % start RISE
%     rise_startup()
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
    % we just give the start date, RISE automatically understand that we
    % are dealing with quarterly data by the format startdate
    mydata.(dataList{id,1})=rise_time_series(startdate,... start date
        tmp.qdatae(:,id+1),... the data
        dataList{id,2});
end

%% do further transformations if necessary

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
        'rise_save_macro',true,... % we ask rise to write the expanded model to disk
        'data',mydata... % we may assign the data now or later
        );
    % a model with multiple files inserted can be difficult to read. The
    % expanded model could be useful for understanding what RISE does and
    % for debugging purposes. The expanded model contains all the details
    % of the individual files (without the comments)
end

%% We estimate the models or filter them directly
profile off
profile on
close all
% if we have the parallel computing toolbox, we can estimate all models in
% one go
for imod=1:nmodels %
    % replace "for" by "parfor" if you want to use parallel computation
    disp('*--------------------------------------------------------------*')
    disp(['*-----Estimation of ',model_names{imod},' model-------*'])
    disp('*--------------------------------------------------------------*')
    estim_models{imod}=estimate(estim_models{imod},'optimizer','fmincon');
    close all
end
profile off
profile viewer

%% Simulation of the posterior distribution
% For each model, the draws are stored in folder
% estim_models{imod}.folders_paths.simulations
mcmc_number_of_simulations=500;
for imod=1:nmodels
    posterior_simulator(estim_models{imod},...
        'mcmc_number_of_simulations',mcmc_number_of_simulations);
end

%% for each set of simulations we compute irfs and plot them
% Since we have a switching model, it would be costly to construct
% generalized impulse responses, for each parameter vector as the
% generalized irfs require a great deal of sampling themselves.
% We take a shortcut for now: we set the number of draws for the
% generalized irfs to 1.
ci=65/100; % confidence region
irf_periods=20; % length of irfs
shock_list=estim_models{1}.exogenous.name;
shock_list_tex=estim_models{1}.exogenous.tex_name;
var_list={'PAI','Y','I'}';
nvar=numel(var_list);
nshocks=numel(shock_list);
locs=locate_variables(var_list,estim_models{1}.endogenous.name);
var_list_tex=estim_models{1}.endogenous.tex_name(locs);
ndraws=mcmc_number_of_simulations;
for imod=1:nmodels
    % go into the folder and load the draws
    W=what(estim_models{imod}.folders_paths.simulations);
    % get the mat files only
    W=W.mat;
    itercount=0;
    % initialize an array for irfs
    this_model=nan(irf_periods+1,nvar,nshocks,ndraws);
    for ifile=1:numel(W)
        tmp=load([estim_models{imod}.folders_paths.simulations,filesep,W{ifile}]);
        for iparam=1:size(tmp.Params,2)
            % get a parameter vector
            params=tmp.Params(:,iparam);
            % push each parameter inside the model and re-solve it
            m=solve(estim_models{imod},'evaluate_params',params);
            itercount=itercount+1;
            if itercount==size(this_model,4)
                % expand the array if necessary
                this_model(:,:,:,end+(1:100))=nan;
            end
            % compute irfs based on the new parameter
            myirfs=irf(estim_models{imod},'irf_type','girf',...
                'girf_draws',1,...
                'irf_periods',irf_periods);
            % extract the irfs
            for ishock=1:nshocks
                this_shock=shock_list{ishock};
                for ivar=1:nvar
                    this_var=var_list{ivar};
                    this_model(:,ivar,ishock,itercount)=...
                        double(myirfs.(this_shock).(this_var));
                end
            end
        end
    end
    
    % sort the irfs, making sure we don't have the nan
    this_model=sort(this_model(:,:,:,1:itercount),4);
    % compute the mean
    this_mean=mean(this_model,4);
    % get the quantiles
    left=round(0.5*(1-ci)*itercount);
    right=round(0.5*(1+ci)*itercount);
    % plot the irfs
    for ishock=1:nshocks
        this_shock=shock_list_tex{ishock};
        figure('name',['irfs for to a ',this_shock,' in model ',model_names{imod}])
        for ivar=1:nvar
            this_var=var_list_tex{ivar};
            subplot(3,1,ivar)
            plot([this_mean(:,ivar,ishock),squeeze(this_model(:,ivar,ishock,[left,right]))])
            title(this_var)
        end
    end
end

