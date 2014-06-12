%% housekeeping
clear all
close all
clc
fclose('all');
%% choice of derivatives

derivatives='symbolic';
%%
addpath('C:\Users\Junior\Documents\GitHub\RISE_toolbox'),rise_startup()

%% read the model
lwz=rise('lwz_model');
%% get the list of the observed variables

obsList={lwz.varobs.name};

% //options_.Harvey_scale_factor = 10; //100;
% 
% estimation(datafile=data_diff_FHFA,mode_check,mode_compute=5,lik_init=2,mh_replic=0,mh_nblocks=2,presample = 3,mh_jscale=0.20,mh_drop=0.2,optim=('Display','iter','MaxFunEvals',1500),smoother,forecast=20) DLogQl  DLogQ  DLogC DLogI DLogB LogL Ql I B C Y;


%% load the data
run data_diff_FHFA
db=struct();
for iobs=1:numel(obsList)
   db.(obsList{iobs})=ts('1975q1',eval(obsList{iobs})); 
end

%% estimate the model
close all,clc
profile off
profile on
switch derivatives
    case 'symbolic'
        lwz=estimate(lwz,'data',ts.collect(db),'optimizer','fmincon')
    case 'analytic'
        lwz=estimate(lwz,'data',ts.collect(db),'optimizer','fmincon',...
            'derivatives','numerical')
    otherwise
end
profile off
profile viewer
