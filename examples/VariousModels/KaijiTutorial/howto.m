%% housekeeping
close all
clear all
clc
%% load the rise paths
addpath('C:\Users\Junior\Documents\GitHub\RISE_toolbox')

%% load RISE
rise_startup();

%% read the model
kaiji=rise('kaiji1207');

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
F=FP/P; % P should always be invertible unless collinearity
% recover the shock impact on the states
K=R(state_cols,:);
%% compute impulse responses
% In general, if you want to see the options for a particular method,
% create an empty rise object. e.g. tmp=rise.empty(0); then you can go
% ahead and call the specific function you want to use on that empty object
% and it will list out the various options and the defaults of those
% options. e.g. irf(tmp)
simple_irfs=irf(kaiji,'irf_periods',20);
%% construct a vector of models with different anticipation options
kaiji_unant=kaiji.set_options('irf_anticipate',false);
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
    figure('name',['orthogonalized shocks to ',shock]);
    for ivar=1:numel(var_list)
        endovar=var_list{ivar};
        tmp=[myirfs{1}.(shock).(endovar),myirfs{2}.(shock).(endovar)];
        subplot(3,3,ivar)
        plot(tmp,'linewidth',2)
        title(endovar)
        if ivar==1
            legend('anticipated','unanticipated')
        end
    end
end
%% estimation

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

kaiji=estimate(kaiji,'data',rise_time_series.collect(mydata))