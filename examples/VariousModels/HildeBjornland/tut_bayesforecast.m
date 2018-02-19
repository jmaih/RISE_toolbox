% housekeeping
close all
clear
clc
%% load estimated models

load tut_bayestimation

load tut_data

endog=models.ve.endogenous;

%% posterior sampling of parameters
params=struct();
params.ve=models.ve.sampler(1000);
params.ve_lr=models.ve_lr.sampler(1000);

%% forecast over all sampled parameters
myfkst=struct();
% date_start=[];
date_start='2003Q1';
myfkst.ve=forecast(models.ve,db,date_start,params.ve);
myfkst.ve_lr=forecast(models.ve_lr,db,date_start,params.ve_lr);

%% set environment

ci=[30,50,68,90];

%% IRFs plots

modelnames=fieldnames(models);

for jj=1:numel(modelnames)
    
    modname=modelnames{jj};
    
    figure('name',['model (',modname,') Projections des variables Norvegiennes']);
    
    for ii=1:numel(endog)
        
        subplot(3,2,ii)
        
        d=myfkst.(modname).(endog{ii});
        
        out=fanchart(d,ci);
        
        plot_fanchart(out)
        
        %     hold on
        %
        %     plot('1995:2004',db.(endog{ii}))
        
        title(tex.(endog{ii}))
        
        axis tight
        
    end
    
    xrotate(45)
    
end


