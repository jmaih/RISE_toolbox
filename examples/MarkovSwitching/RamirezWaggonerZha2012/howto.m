%% Housekeeping
clear
close all
clc

%% read the models and their calibrations

m0=rise('frwz_nk');

%% solve the model

m=solve(m0,'steady_state_unique',true);

part_list={};

m2=solve(m0,'solve_perturbation_type',{'frwz',part_list});


%% print results
m.print_solution()

%% print solution for a subset of variables only
m.print_solution({'PAI','Y'})

%% Alternative calibrations

cal_2=struct('psi_a_2', 0.7);

m2=set(m,'parameters',cal_2);

%% construct a vector of models

bigm=[m,m2];

%% compute impulse responses for all models simultaneously

myirfs=irf(bigm,'irf_periods',20);

%% plot the impulse responses
close all
var_list=m.endogenous.name;
figure('name','Impulse responses to a monetary policy shock');
for ii=1:numel(var_list)
    subplot(3,1,ii)
    reg1=myirfs.EPS_R.regime_1.(var_list{ii});
    reg2=myirfs.EPS_R.regime_2.(var_list{ii});
    plot([reg1,reg2]);
    title(var_list{ii})
    if ii==1
        legend({'reg1-m1','reg1-m2','reg2-m1','reg2-m2'})
    end
    axis tight
end

