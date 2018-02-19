% housekeeping
close all
clear
clc
%% load estimated models

load tut13_bayestimation

load tut01_data

endog=models.ve.endogenous;

%% posterior sampling of parameters
params=struct();
params.ve=models.ve.estim_.sampler(1000);
params.ve_lr=models.ve_lr.estim_.sampler(1000);

%% forecast over all sampled parameters
clc
myfkst=struct();
% date_start=[];
date_start='2003Q1';
nsteps=12;
shock_uncertainty=false;
Rfunc=[];

conditions=struct();
conditions.rtwi={'2003Q1','2004Q4'}; % range over which we want to condition

myfkst.ve=forecast(models.ve,db,date_start,params.ve,nsteps,...
    shock_uncertainty,Rfunc,conditions);

myfkst.ve_lr=forecast(models.ve_lr,db,date_start,params.ve_lr,nsteps,...
    shock_uncertainty,Rfunc,conditions);

%% set environment

ci=[30,50,68,90];

%% Forecast plots

modelnames=fieldnames(models);

for jj=1:numel(modelnames)
    
    modname=modelnames{jj};
    
    figure('name',['model (',modname,') Forecasts of Norwegian Data']);
    
    for ii=1:numel(endog)
        
        subplot(3,2,ii)
        
        d=myfkst.(modname).(endog{ii});
        
        out=fanchart(d,ci);
        
        plot_fanchart(out,[244, 122, 66]/255)%'c'
        
        %     hold on
        %
        %     plot('1995:2004',db.(endog{ii}))
        
        title(tex.(endog{ii}))
        
        axis tight
        
    end
    
    xrotate(45)
    
end
