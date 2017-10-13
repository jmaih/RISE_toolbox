%% housekeeping 
clear
clc
close all
%% rise the model

m=rise('ig2017');

%% assign parameters and steady state file

m=set(m,'parameters',paramfile(),'steady_state_file','ssfile');

%% solve the occbin model
% assign each restriction to a markov chain
restr_map=struct('debt',1,'zlb',2);

m=solve(m,'solve_occbin',{1,restr_map,m.markov_chains.regimes},...
    'steady_state_unique',true,... % impose a unique steady state
    'steady_state_imposed',true);

%% simulate the model (honoring the constraints)
rng(1900)
mysims=simulate(m,'simul_honor_constraints',true,'simul_periods',300);

%% plot the simulations
close all
figure();
vnames={'b','maxlev','debt','lm','zlb','r','rnot'};
for ii=1:numel(vnames)
    v=vnames{ii};
    subplot(3,3,ii)
    plot(mysims.(v),'linewidth',2)
%     plot('115:200',mysims.(v),'linewidth',2)
    title(v)
    axis tight
end