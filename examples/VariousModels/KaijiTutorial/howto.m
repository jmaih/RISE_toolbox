%% housekeeping
close all
clearvars
clc

%% read the model
m=rise('kaiji1207_exp');

%% alternatively, solve the model directly

m=solve(m);

%% print the solution of the model

m.print_solution
% alternative call
% print_solution(kaiji)

%% print solution of a subset of variables

m.print_solution({'C','H','K'})

%% doing things the RBC way: waste of time
% RISE write the solution of the whole system as: X_t=T*X_{t-1}+R*e_t. The
% RBC guys write solutions as S_t=P*S_{t-1}+K*e_t and Y_t=F*S_t . Below, we
% show how to recover matrices P, K and F.
% collect the T and R matrices
[T,R]=load_solution(m,'iov');
T=T{1}; R=R{1};
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

allvars=m.endogenous.name;

state_variables=allvars(state_cols);

control_variables=allvars(control_cols);

%% compute impulse responses
% In general, if you want to see the options for a particular method,
% create an empty rise object. e.g. tmp=rise.empty(0); then you can go
% ahead and call the specific function you want to use on that empty object
% and it will list out the various options and the defaults of those
% options. e.g. irf(tmp)
simple_irfs=irf(m,'irf_periods',20);
%% construct a vector of models with different anticipation options
kaiji_unant=m.set('irf_anticipate',false);
myvector=[m,kaiji_unant];

%% change some option in the vector
% here we set:
% 1- the number of irf periods. If we don't the default of 40 will be used
% 2- we assume that we know as of today that shocks will hit 3 periods from
% now. And we choose to act on that information (anticipated) or disregard
% the information (unanticipated)
% 3- we can also choose the sign of the shocks. By default, the shocks will
% be positive.
myirfs=irf(myvector,'irf_periods',20,'solve_shock_horizon',3,'irf_shock_sign',1);

%% plot the irfs
shock_list=m.exogenous.name;
tex=get(m,'tex');
% var_list=allvars;
% just re-ordering the variables in a way that I like for the plotting
var_list={'A','B','D','PSI','C','K','H','R'}; 
for ishock=1:numel(shock_list)
    shock=shock_list{ishock};
    figure('name',['IRFs to a ',tex.(shock), 'shock']);
    for ivar=1:numel(var_list)
        endovar=var_list{ivar};
        subplot(3,3,ivar)
        plot(myirfs.(shock).(endovar),'linewidth',2)
        title(tex.(endovar))
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
    mydata.(names{ii})=ts(start_date,datta(:,ii+1));
end

%% plot the data
figure('name','Observed data');
for ii=1:numel(names)
    subplot(3,1,ii)
    plot(mydata.(names{ii}))
    title(names{ii});
end

%% estimate the model

m=estimate(m,'data',mydata);

%% historical decomposition of shocks
histdec=historical_decomposition(m);

%% plot the decomposition
figure('name','historical decomposition of shocks and initial conditions')
for ivar=1:numel(var_list)
    vname=var_list{ivar};
    subplot(3,3,ivar)
    plot_decomp(histdec.(vname))
    title(tex.(vname))
    if ivar==1
        contrib_names=histdec.(vname).varnames;
        for jj=1:numel(contrib_names)
            if ~isfield(tex,contrib_names{jj})
                continue
            end
            contrib_names{jj}=tex.(contrib_names{jj});
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
    %,{'EPS_PSI','EPS_B','EPS_D'}
    [counterf,actual]=counterfactual(m,[],1,shock_list{ishock});
    figure('name',['Counterfactual: ',tex.(shock_list{ishock}),' shock only'])
    for ivar=1:numel(var_list)
        vname=var_list{ivar};
        subplot(3,3,ivar)
        plot([actual.(vname),counterf.(vname)])
        title(tex.(vname))
        if ivar==1
            legend({'actual','counterfactual'})
        end
    end
end

%% counterfactual for a subset of shocks
close all
[counterf,actual]=counterfactual(m,[],1,{'EPS_PSI','EPS_B','EPS_D'});
figure('name','Counterfactual: EPS_PSI and EPS_B shock only');
for ivar=1:numel(var_list)
    vname=var_list{ivar};
    subplot(3,3,ivar)
    plot([actual.(vname),counterf.(vname)])
    title(tex.(vname))
    if ivar==1
        legend({'actual','counterfactual'})
    end
end

%% variance decomposition
vardec=variance_decomposition(m);

%% plot the decomposition
close all
figure('name','Variance decomposition of shocks')
for ivar=1:numel(var_list)
    vname=var_list{ivar};
    subplot(3,3,ivar)
    plot_decomp('0:50',vardec.conditional.(vname))
    title(tex.(vname))
    if ivar==1
        contrib_names=vardec.conditional.(vname).varnames;
        for jj=1:numel(contrib_names)
            if ~isfield(tex,contrib_names{jj})
                continue
            end
            contrib_names{jj}=tex.(contrib_names{jj});
        end
        hleg=legend(contrib_names,...
            'Location','BestOutside','orientation','horizontal');
        pp=get(hleg,'position');
        pp(1:2)=0;
        set(hleg,'position',pp)
    end
end
