%% housekeeping
close all
clear all
clc
%% add the necessary paths
setpaths

%% get the data
[datt,navn]=xlsread('DataNEMO');
startdate=[int2str(datt(1,1)),'Q',int2str(datt(1,2))];
rawdata=rise_time_series(startdate,datt(:,3:end),char(navn(3:end)));
obs_names=char('DPQ_P_NW','D_GDP_NW','RN3M_NW');
cond_obs_names=char('RN3M_NW');
data=window(rawdata,startdate,'',obs_names);

%% Two-Markov-chain model object
TwoChain=rise('Canonical_2_MarkovChains','data',data,...
    'optimset',optimset('MaxNodes',20,'TolFun',sqrt(eps),'MaxTime',2*60*60,...
    'MaxIter',2000),'check_stability',false,'solver','functional_iteration');%newton_kronecker_iteration
%%
TwoChain=TwoChain.estimate('optimizer','bee_gate');%
% TwoChain=TwoChain.estimate(1);%
% check_optimum(TwoChain)
