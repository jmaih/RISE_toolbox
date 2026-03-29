%% housekeeping
close all
clear
clc

%% Read the models

m1=rise('usmodel_b','rise_flags',{'policy','taylor'});
m2=rise('usmodel_b','rise_flags',{'policy','optimal'});
m2=set(m2,'fix_point_maxiter',5000,... % jack up the number of iterations to increase the probability of solving
    'solve_check_stability',false);% save some time by avoiding the checking of stability all the time

%% load the data and construct the time series
dat=load('usmodel_data');
vnames=fieldnames(dat);
data=[];
for jj=1:numel(vnames)
    data=[data,dat.(vnames{jj})];
end
data_start='1947q3';
data=ts(data_start,data,vnames);


%% estimating the models

nw=utils.parallel.get_number_of_workers();

if nw==0
    mypool=parpool(2);
end

sw=estimate([m1,m2],'data',data,...
    'estim_start_date',obs2date(data_start,71),... % retrieve the date of the 71st observation
    'kf_presample',4); % <--- sw=sw.estimate;

if nw==0
    delete(mypool)
end