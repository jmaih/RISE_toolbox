%% housekeeping
close all
clear
clc
%% load the toolbox
rise_startup()
%% load the data and construct the time series
dat=load('usmodel_data');
vnames=fieldnames(dat);
data=[];
for jj=1:numel(vnames)
    data=[data,dat.(vnames{jj})];
end
data_start='1947q3';
data=rise_time_series(data_start,data,vnames);
%% read the model and assign the data
hist_start=rise_date(data_start);
sw=rise('usmodel_switch_ssfile','data',data,...
    'steady_state_file','usmodel_lc_steady_state',...
    'estim_start_date',hist_start.observation_2_date(71),...
    'kf_presample',4,...
    'lc_MaxIter',5000,... % jack up the number of iterations to increase the probability of solving
    'check_stability',false); % save some time by avoiding the checking of stability all the time
    

%% estimating the model
profile on
sw=sw.estimate;
profile off
profile viewer