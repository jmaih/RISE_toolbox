%% housekeeping
close all
clear
clc
%% load the data

load tut01_data

load tut02_estimation

load tut04_identification

%% pick a model

% play with different models to see interesting results!!!

model='ve_lr'; % ve_lr

endog=models.(model).endogenous;

%% choose identification

Rfunc=Rfunc3;

%% compute decompositions

params=[];

vd=variance_decomposition(models.(model),params,Rfunc3.(model),shock_names);

%% plot decompositions

range='0:50'; % pick a range for the plots

figure('name','Variance Decomposition');

for iv=1:numel(endog)
    
    d=vd.conditional.(endog{iv});
    
    subplot(3,2,iv)
    
    plot_decomp(range,d)
    
    if iv==1
        
        legend(shock_names,'location','SE',...
            'Orientation','horizontal')
        
    end
    
    title(tex.(endog{iv}))
    
end