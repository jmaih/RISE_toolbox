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
data=ts(data_start,data,vnames);
%% read the model and assign the data
sw=rise('usmodel_switch_ssfile','data',data,...
    'steady_state_file','usmodel_lc_steady_state',...
    'estim_start_date',obs2date(data_start,71),... % retrieve the date of the 71st observation
    'kf_presample',4,...
    'lc_MaxIter',5000,... % jack up the number of iterations to increase the probability of solving
    'solve_check_stability',false); % save some time by avoiding the checking of stability all the time
    

%% estimating the model
sw=estimate(sw); % <--- sw=sw.estimate;