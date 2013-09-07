%% housekeeping
close all
clear
clc
%%
cd C:\Users\Junior\Dropbox\RISE\RISE_Toolbox\examples\VariousModels\SmetsWouters
%% add the necessary paths
addpath C:\Users\Junior\Dropbox\RISE\RISE_Toolbox
rise_startup();
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
sw=rise('usmodel','data',data,...
    'estim_start_date',hist_start.observation_2_date(71),...
    'presample',4,...
    'steady_state_file','usmodel_steadystate',...
    'solver','msre_aim',...
    'derivatives','numerical');

%% estimating the model
profile on
sw=sw.estimate;
profile off
profile viewer