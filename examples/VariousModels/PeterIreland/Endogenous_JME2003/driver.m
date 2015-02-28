%% housekeeping
clear all
close all
clc
%% RISE the model
m=rise('ireland2003');

%% load the data
mydat=load('cimpr.dat');
bigt = size(mydat,1);

ct=mydat(:,1);
it=mydat(:,1);
mt=mydat(:,1);
trend = (1:bigt)';
xxx = [ ones(bigt,1),trend];

betac = xxx\ct;
betai = xxx\it;
betam = xxx\mt;

ct = ct - betac(2)*trend;
it = it - betai(2)*trend;
mt = mt - betam(2)*trend;

mydat(:,1:3)=[ct,it,mt];
data=ts('1959Q1',mydat,{'LC','LI','LM','LPI','LR'});
data=pages2struct(data);
data.TREND=ts(data.LC.start,trend);
%%
ms=estimate(m,'data',data,'kf_tol',0);
%% Redo filtering?

mfilt=filter(ms,'data',data,'kf_tol',0);