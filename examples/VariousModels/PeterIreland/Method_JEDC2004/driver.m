%% housekeeping
clear all
close all
clc
%% RISE the model
m=rise('ireland2004');

%% load the data
mydat=load('ych.dat');
data=ts('1948Q1',log(mydat),{'LY','LC','LH'});
data=pages2struct(data);
data.TREND=ts(data.LY.start,(1:data.LY.NumberOfObservations).');
%%
ms=estimate(m,'data',data,'kf_tol',0);
%% Redo filtering?

mfilt=filter(ms,'data',data,'kf_tol',0);