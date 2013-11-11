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
sw=rise('usmodel','data',data,...
    'estim_start_date',obs2date(data_start,71),... % retrieve the date of the 71st observation
    'kf_presample',4,...
    'steady_state_file','usmodel_steadystate');

%% estimating the model
sw=estimate(sw); % <--- sw=sw.estimate;
