%% housekeeping
clear classes 
close all 
clc
%% load RISE
do_set_paths=true;
if do_set_paths
    addpath C:\Users\Junior\Documents\GitHub\RISE_toolbox
    rise_startup()
end
%% Load and transform data
db=mertens_ravn_data();
%% add a trend, a quadratic trend and a dummy for 1975q2 
% this is econometric nonsense
vnames=fieldnames(db);
prototype=db.(vnames{1});
nobs=prototype.NumberOfObservations;
start=prototype.start;
finish=prototype.finish;
mytrend=(1:nobs)';
% I do not like huge numbers. So I recast everything in the interval [0,1].
% the drawback is that this is not valid out of sample. If you do, then
% comment the line below. In any case, pick your poison :)
mytrend=(mytrend-min(mytrend))/(max(mytrend)-min(mytrend)); 
db.trend=rise_time_series(prototype.TimeInfo,mytrend);
db.qtrend=rise_time_series(prototype.TimeInfo,mytrend.^2);
db.dum1975q2=rise_time_series.dummy(prototype.TimeInfo,[],'1975q2');
%% create a svar as a rise object
mysvar=struct('model','svar',...
    'var',{{'RGDP',['$',db.RGDP.varnames{1},'$'],...
    'GOV',['$',db.GOV.varnames{1},'$'],...
    'APITR',['$',db.APITR.varnames{1},'$']}},...
    'varexo',{{'trend','qtrend','dum1975q2'}});
% push the order of the VAR directly, and override the default of 4. You
% could always do it later as well
obj=rise(mysvar,'data',db,'svar_lags',4,'svar_restarts',10);%
%% push the restrictions into the object
% in this case the restrictions are written in an m-file. Without
% arguments, the function should return the list of the endogenous
% variables as ordered in the impact matrix. With one argument, the
% parameters to estimate, the function should return R, the impact matrix.
% In all cases, the function should always return in addition to R, the
% lower and the upper bounds for the parameters to estimate.
obj=set_options(obj,'svar_restrictions',@blanchard_perotti_restrictions);
%% Solve all 4 models simultaneously
clc
objs=solve(obj,'debug',true);
%% compute irfs for all 3 models simultaneously
testirfs=irf(obj,'irf_periods',20,'irf_shock_sign',-1);
%% plot the irfs
figure('name','impulse response functions')
endo_names=fieldnames(testirfs.(obj(1).varexo(1).name));
varexo_det_names={obj(1).varobs_exo.name};
iter=0;
for iname=1:numel(endo_names)
    v=endo_names{iname};
    vloc=find(strcmp(v,{obj(1).varendo.name}));
    vtex=obj(1).varendo(vloc).tex_name;
    for ishock=1:sum(obj(1).exogenous.number)
        shock_name=obj(1).varexo(ishock).name;
        if ismember(shock_name,varexo_det_names)
            continue
        end
        shock_texname=obj(1).varexo(ishock).tex_name;
        iter=iter+1;
        subplot(3,3,iter)
        plot(testirfs.(shock_name).(v),'linewidth',2)
        if ishock==1
            ylabel(vtex)
        end
        if iname==1
            title(shock_texname,'interpreter','none')
        end
        if iname==1 && ishock==1
            leg=legend({obj.filename});
        end
    end
end

