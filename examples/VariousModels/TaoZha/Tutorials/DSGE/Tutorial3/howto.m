%%
% RISE Tutorial by Dr. Tao Zha
%% housekeeping
clear all
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
        'data',mydata... % we may assign the data now or later
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
for imod=1:nmodels %
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
mcmc_number_of_simulations=500;
do_store_mcmc_to_disk=false;
postSims=cell(1,nmodels);
obj=cell(1,nmodels);
for imod=1:nmodels
    % replace "for" by "parfor" if you want to use parallel computation
    disp('*--------------------------------------------------------------*')
    disp(['*----------- MCMC for ',model_names{imod},' model-------------*'])
    disp('*--------------------------------------------------------------*')
    if do_store_mcmc_to_disk
        % In this case, the simulations are either stored on disk. You may
        % want to do this if you need a very large number of draws. 
        %-----------------------------------------------------------------
        obj{imod}=posterior_simulator(estim_models{imod},'mcmc_number_of_simulations',...
            mcmc_number_of_simulations);
    else
        % Here the simulations are stored in variable, which here we call
        % postSims.
        %-----------------------------------------------------------------
        [obj{imod},postSims{imod}]=posterior_simulator(estim_models{imod},'mcmc_number_of_simulations',...
            mcmc_number_of_simulations);
    end
end
%% Marginal Likelihood
% choose the algorithm for computing the marginal data density in case you
% want to compute it. The algorithms for the computation of MDD are still
% work in progress
%-------------------------------------------------------------------------
lmdd_algo={'mhm','chib_jeliazkov'};
choice=lmdd_algo{1};
lmdd=cell(1,nmodels);
do_log_marginal_data_density=false;
if do_log_marginal_data_density
    for imod=1:nmodels
    % replace "for" by "parfor" if you want to use parallel computation
    disp('*--------------------------------------------------------------*')
    disp(['*--- Marginal data density for ',model_names{imod},' model---*'])
    disp('*--------------------------------------------------------------*')
        if do_store_mcmc_to_disk
            lmdd{imod}=log_marginal_data_density(obj{imod},choice);
        else
            lmdd{imod}=log_marginal_data_density(obj{imod},choice,postSims{imod});
        end
    end
end

%% posterior distribution of parameters
% to be added later

%% draw random parameters from the simulation above and compute irfs
% choose the number of parameter draws for the irfs or a similar exercise
%------------------------------------------------------------------------
myreplications=200;

batch_irfs=cell(1,nmodels);
parfor imod=1:nmodels
    % replace "for" by "parfor" if you want to use parallel computation
    disp('*--------------------------------------------------------------*')
    disp(['*-------- Bayesian IRFs for ',model_names{imod},' model--------*'])
    disp('*--------------------------------------------------------------*')
    ireplic=0;
    while ireplic<myreplications
        if do_store_mcmc_to_disk
            [draw,obj_with_draw]=draw_parameter(obj{imod});
        else
            [draw,obj_with_draw]=draw_parameter(obj{imod},postSims{imod});
        end
        % solve the model
        %----------------
        [obj_with_draw,retcode]=solve(obj_with_draw);
        % if there is no problem proceed. else draw a new parameter. But
        % since the parameters have been mcmc-ed. No problem is to be
        % expected. So this is done for comprehensiveness since the
        % draw_parameter function also allows you to draw directly from the
        % prior distribution too 
        %------------------------------------------------------------------
        if retcode
            continue
        end
        ireplic=ireplic+1;
        % now you can compute your irfs based on that single draw
        %--------------------------------------------------------
        myirfs=irf(obj_with_draw);
        if ireplic==1
            shock_names=fieldnames(myirfs);
            vnames=fieldnames(myirfs.(shock_names{1}));
            regime_names=myirfs.(shock_names{1}).(vnames{1}).varnames;
            for ishock=1:numel(shock_names)
                for ivar=1:numel(vnames)
                    batch_irfs{imod}.(shock_names{ishock}).(vnames{ivar})=[];
                end
            end
        end
        for ishock=1:numel(shock_names)
            for ivar=1:numel(vnames)
                batch_irfs{imod}.(shock_names{ishock}).(vnames{ivar})=cat(3,...
                    batch_irfs{imod}.(shock_names{ishock}).(vnames{ivar}),...
                    double(myirfs.(shock_names{ishock}).(vnames{ivar})));
                if ireplic==myreplications
                    % put back into time series format
                    %---------------------------------
                    batch_irfs{imod}.(shock_names{ishock}).(vnames{ivar})=...
                        ts(0,...
                        batch_irfs{imod}.(shock_names{ishock}).(vnames{ivar}),...
                        regime_names);
                end
            end
        end
    end
end
%% prepare a report but now instead of publishing directly, we will add the irfs below
close all
rep_data=struct();
myfields={'name','title','address','author','email','abstract'};

rep_data.name='RudebuschSvensson'; % Name under which the report will be saved
rep_data.title='Switches in the US Macroeconomic Data using the Rudebusch-Svensson VAR model';
rep_data.address='Your institution'; 
rep_data.author='first\_name last\_name'; % your name 
rep_data.email='first\_name@last\_name.com'; % your email
myabstract={['This report investigates switches in the parameters of a ',...
    ' simple VAR model by Rudebusch and Svensson (1999) estimated on US data. 4 variants of the model',...
    'are estimated: (i) the first model has constant parameters; ',...
    '(ii) the second model allows for switches in the policy parameters only; ',...
    '(iii) the third specification allows for switches in volatility only ;',...
    '(iv) the fourth variant allows for independent switches in both parameters and the ',...
    'volatility of shocks.']
    ' ' % leave a blank space to start the next sentence in a new line
    'We find ample evidence in favor of switching parameters... '
    };
rep_data.abstract=myabstract;

xrep=Tao_report([obj{:}],rep_data);
%% distribution of the irfs
% set the confidence regions for your irfs
%-----------------------------------------
confi_irfs=[.3,.5,.7,.9];
% choose the colors for the fan charts
%-------------------------------------
mycolors={'nb','y','g','b','m','r'};
mycolors=mycolors{1};
% plot a shorter horizon for the irfs if you want. 40 is the default in
% RISE
%-----------------------------------------------------------------------
irf_plot_horizon=40;
myrange=sprintf('0:%0.0f',irf_plot_horizon);
for imod=1:nmodels
    disp('*--------------------------------------------------------------*')
    disp(['*-------- Plotting IRFs for ',model_names{imod},' model--------*'])
    disp('*--------------------------------------------------------------*')
    shock_names=fieldnames(batch_irfs{imod});
    % Find the list of the endogenous variables that are not auxiliary
    % (auxiliary variables are variables created so that the model has a
    % maximum number of lags of 1 and similarly for the leads if the model
    % is forward looking
    %-------------------------------------------------------------------
    vnames=get(obj{imod},'endo_list(~auxiliary)'); 
    locs=locate_variables(vnames,obj{imod}.endogenous.name);
    vtex_names=obj{imod}.endogenous.tex_name(locs);
    nvars=numel(vnames);
    regime_names=batch_irfs{imod}.(shock_names{1}).(vnames{1}).varnames;
    for ishock=1:numel(shock_names)
        figh=nan(1,numel(regime_names));
        for iregime=1:numel(regime_names)
            figh(iregime)=figure('name',[model_names{imod},':: Posterior ',...
                'Impulse responses to a ',shock_names{ishock},' shock in ',...
                regime_names{iregime}]);
        end
        for ivar=1:numel(vnames)
            mydata=batch_irfs{imod}.(shock_names{ishock}).(vnames{ivar});
            for iregime=1:numel(regime_names)
                out=fanchart(mydata(myrange,mydata.varnames{iregime},:),confi_irfs);
                figure(figh(iregime))
                % you may want to change the line below if you want a
                % different configuration for your plots
                %----------------------------------------------------------
                subplot(nvars,1,ivar)
                h=plot_fanchart(out,mycolors);
                title(vtex_names{ivar})
                axis tight
                if ivar==numel(vnames)
                    legend(char(num2str(flipud(confi_irfs(:))),'mean'),'location','SW','orientation','horizontal')
                end
            end
        end
        for ifig=1:numel(figh)
            xrep.figure('name',figh(ifig),'title',get(figh(ifig),'Name'))
            xrep.pagebreak();
        end
    end
end

%% now we publish the report
publish(xrep)
close all


