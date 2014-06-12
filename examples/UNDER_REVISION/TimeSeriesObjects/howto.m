%% create empty time series

t0a=ts.empty(0,1);
 % or simply
t0b=ts;
%% Generate a 6 random matrix which will be our dataset
% we create time series with 20 observations
smpl=20;

RM=rand(smpl,6);
% we want to turn those data into time series
%% create a time series at different frequencies, with starting date 1990
% for simplicity, assume 
%  column 1 of RM represents yearly data of a variable whose name is v1
%  column 2 of RM represents bi-annual data of a variable whose name is v2
%  column 3 of RM represents quarterly data of a variable whose name is v3
%  column 4 of RM represents monthly data of a variable whose name is v4
%  column 5 of RM represents weekly data of a variable whose name is v5
%  column 6 of RM represents daily data of a variable whose name is v6

% at yearly frequency: 2 options
v1a=ts(1990,RM(:,1));
v1b=ts('1990',RM(:,1));

% at bi-annual frequency
v2=ts('1990H1',RM(:,2));

% at quarterly frequency
v3=ts('1990Q1',RM(:,3));

% at monthly frequency
v4=ts('1990M1',RM(:,4));

% at weekly frequency: 2 options
v5a=ts('1990M1D5W',RM(:,5));
v5b=ts('19900105W',RM(:,5));

% at daily frequency: 2 options
v6a=ts('1990M1D1',RM(:,6));
v6b=ts('19900101D',RM(:,6));
%% display a time series
% let's display the last one
v6b
% let's display the data
v6b.data
%% create leads and lags of say the annual series
v1a_lead=v1a(+4);
v1a_lag=v1a(-3);
%% collect all those annual data in one time series
v1c=[v1a,v1a_lag,v1a_lead];
%% display those data
v1c.data
%% plot them
figure;
plot(v1c)
figure;
bar(v1c)
%% create a database of the 3 series with names
% get the data as a matrix
annual_data=double(v1c);
% the vector of dates is already computed for us
v1d=ts(v1c.TimeInfo,annual_data,{'current','lags','leads'})
% plot again and add a legend
figure;
plot(v1d)
legend({'current','lags','leads'})
%% let's do some arithmetics
% illustrating the use of overloaded operations + - / *
v1_arit=v1a-v1a(-1)+v1a*log(v1a(+2))/exp(v1a);
%% check the normality of our annual series
jbtest(v1a)
%% skewness and kurtosis
skewness(v1a)
kurtosis(v1a)
%% estimate an ar(3) process: a constant is included at the end by default
v1a.ar(3)
%% regression form
regress(v1a,v1a(-1)&v1a(-2)&v1a(-3))
% add a vector of ones at the end to get the same results as the ar(3)
% above
regress(v1a,ones(v1a(-1)&v1a(-2)&v1a(-3)))
%% collect time series of different frequencies
% this operation will convert all the series to the lowest frequency and
% then put them in the same database
collection=ts.collect({'v1',v1a},{'v2',v2},{'v3',v3},{'v4',v4},{'v5',v5a},{'v6',v6a});


