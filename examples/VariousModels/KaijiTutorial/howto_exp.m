%% housekeeping
close all
clear all
clc
%% load the rise paths
addpath('C:\Users\Junior\Documents\GitHub\RISE_toolbox')

%% load RISE
rise_startup();

%% read the model
kaiji=rise('kaiji1207_exp');

%% get the first-order perturbation (approximation) of the model

kaiji=evaluate(kaiji);

%% alternatively, solve the model directly

kaiji=solve(kaiji);

%% print the solution of the model

kaiji.print_solution
% alternative call
% print_solution(kaiji)

%% print solution of a subset of variables

kaiji.print_solution({'C','H','K'})

%% doing things the RBC way: waste of time
% RISE write the solution of the whole system as: X_t=T*X_{t-1}+R*e_t. The
% RBC guys write solutions as S_t=P*S_{t-1}+K*e_t and Y_t=F*S_t . Below, we
% show how to recover matrices P, K and F.
% collect the T and R matrices
T=kaiji.T;
R=kaiji.R;
% set some elements to 0 sharp
T(abs(T)<1e-10)=0;
R(abs(R)<1e-10)=0;
% get the state columns
state_cols=any(T); 
% Kaiji wastes time double-checking check that find(any(T)) gives the
% columns of the states 
% get the control columns as the columns that are not states
control_cols=~state_cols;
% now build the P matrix expressing the states as a function of the states
P=T(state_cols,state_cols);
% recover the solution of the controls as a function of states
FP=T(control_cols,state_cols);
F=FP/P; % F=FP*inv(P) P should always be invertible unless collinearity
% recover the shock impact on the states
K=R(state_cols,:);

allvars={kaiji.varendo.name};

state_variables=allvars(state_cols);

control_variables=allvars(control_cols);

%% compute impulse responses
% In general, if you want to see the options for a particular method,
% create an empty rise object. e.g. tmp=rise.empty(0); then you can go
% ahead and call the specific function you want to use on that empty object
% and it will list out the various options and the defaults of those
% options. e.g. irf(tmp)
simple_irfs=irf(kaiji,'irf_periods',20);
%% construct a vector of models with different anticipation options
kaiji_unant=kaiji.set('irf_anticipate',false);
myvector=[kaiji,kaiji_unant];

%% change some option in the vector
% here we set:
% 1- the number of irf periods. If we don't the default of 40 will be used
% 2- we assume that we know as of today that shocks will hit 3 periods from
% now. And we choose to act on that information (anticipated) or disregard
% the information (unanticipated)
% 3- we can also choose the sign of the shocks. By default, the shocks will
% be positive.
myirfs=irf(myvector,'irf_periods',20,'irf_horizon',3,'irf_shock_sign',1);

%% plot the irfs
shock_list={kaiji.varexo.name};
var_list={kaiji.varendo.name};
% just re-ordering the variables in a way that I like for the plotting
var_list={'A','B','D','PSI','C','K','H','R'}; 
for ishock=1:numel(shock_list)
    shock=shock_list{ishock};
    loc=locate_variables(shock,{kaiji.varexo.name});
    figure('name',['IRFs to a ',kaiji.varexo(loc).tex_name, 'shock']);
    for ivar=1:numel(var_list)
        endovar=var_list{ivar};
        subplot(3,3,ivar)
        plot(myirfs.(shock).(endovar),'linewidth',2)
        loc=locate_variables(endovar,{kaiji.varendo.name});
        title({kaiji.varendo(loc).tex_name})
        if ivar==1
            legend('anticipated','unanticipated')
        end
    end
end
%% estimation. Now the DSGE way

%% read the data and create the time series
[datta,names]=xlsread('data.xlsx');
names=names(1,2:end);
start_date='1960';

mydata=struct();
for ii=1:numel(names)
    mydata.(names{ii})=rise_time_series(start_date,datta(:,ii+1));
end

%% plot the data
figure('name','Observed data');
for ii=1:numel(names)
    subplot(3,1,ii)
    plot(mydata.(names{ii}))
    title(names{ii});
end

%% estimate the model

kaiji=estimate(kaiji,'data',rise_time_series.collect(mydata));

%% historical decomposition of shocks
histdec=historical_decomposition(kaiji);

%% plot the decomposition
figure('name','historical decomposition of shocks and initial conditions')
for ivar=1:numel(var_list)
    vname=var_list{ivar};
    subplot(3,3,ivar)
    plot_decomp(histdec.(vname))
    loc=locate_variables(vname,{kaiji.varendo.name});
    title(kaiji.varendo(loc).tex_name)
    if ivar==1
        contrib_names=histdec.(vname).varnames;
        shock_texnames={kaiji.varexo.tex_name};
        locs=locate_variables(contrib_names,{kaiji.varexo.name},true);
        for jj=1:numel(locs)
            if isnan(locs(jj))
                continue
            end
            contrib_names{jj}=kaiji.varexo(locs(jj)).tex_name;
        end
        hleg=legend(contrib_names,...
            'Location','BestOutside','orientation','horizontal');
        pp=get(hleg,'position');
        pp(1:2)=0;
        set(hleg,'position',pp)
    end
end

%% counterfactual: what if only one shock had been alive?
for ishock=1:numel(shock_list)
    [counterf,actual]=counterfactual(kaiji,'counterfact_shocks_db',...
        shock_list{ishock});%,{'EPS_PSI','EPS_B','EPS_D'}
    loc=locate_variables(shock_list{ishock},{kaiji.varexo.name},true);
    figure('name',['Counterfactual: ',kaiji.varexo(loc).tex_name,' shock only'])
    for ivar=1:numel(var_list)
        vname=var_list{ivar};
        subplot(3,3,ivar)
        plot([actual.(vname),counterf.(vname)])
        loc=locate_variables(vname,{kaiji.varendo.name});
        title(kaiji.varendo(loc).tex_name)
        if ivar==1
            legend({'actual','counterfactual'})
        end
    end
end

%% counterfactual for a subset of shocks
[counterf,actual]=counterfactual(kaiji,'counterfact_shocks_db',...
    {'EPS_PSI','EPS_B'});
figure('name','Counterfactual: EPS_PSI and EPS_B shock only');
for ivar=1:numel(var_list)
    vname=var_list{ivar};
    subplot(3,3,ivar)
    plot([actual.(vname),counterf.(vname)])
    loc=locate_variables(vname,{kaiji.varendo.name});
    title(kaiji.varendo(loc).tex_name)
    if ivar==1
        legend({'actual','counterfactual'})
    end
end

%% variance decomposition
vardec=variance_decomposition(kaiji);

%% plot the decomposition
figure('name','Variance decomposition of shocks')
for ivar=1:numel(var_list)
    vname=var_list{ivar};
    subplot(3,3,ivar)
    plot_decomp(vardec.conditional.(vname))
    loc=locate_variables(vname,{kaiji.varendo.name});
    title(kaiji.varendo(loc).tex_name)
    if ivar==1
        contrib_names=vardec.conditional.(vname).varnames;
        shock_texnames={kaiji.varexo.tex_name};
        locs=locate_variables(contrib_names,{kaiji.varexo.name},true);
        for jj=1:numel(locs)
            if isnan(locs(jj))
                continue
            end
            contrib_names{jj}=kaiji.varexo(locs(jj)).tex_name;
        end
        hleg=legend(contrib_names,...
            'Location','BestOutside','orientation','horizontal');
        pp=get(hleg,'position');
        pp(1:2)=0;
        set(hleg,'position',pp)
    end
end
