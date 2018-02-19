%% housekeeping
close all
clear
clc
%% load the data

load tut01_data

load tut02_estimation

load tut09_bootstrap

load tut04_identification

endog=models.ve.endogenous;

%% set environment

ci=[30,50,68,90];

% choose identification
quick=true;

Rfunc_star=Rfunc1;

%% irfs 

model_types=fieldnames(models);

myirfs=struct();

for imod=1:numel(model_types)
    
    if quick
        
        Rfunc=[];
        
    else
        
        Rfunc=Rfunc_star.(model_types{imod});
        
    end
    
    myirfs.(model_types{imod})=irf(models.(model_types{imod}),shock_names,...
        40,params.(model_types{imod}),Rfunc);
    
end

%% IRFs plots

for imod=1:numel(model_types)
    
    the_type=model_types{imod};
    
    for ishock=1:numel(shock_names)
        
        shock=shock_names{ishock};
        
        titel=['Model (',the_type,'): IRFs to a ',shock,' shock'];
        
        figure('name',titel);
        
        for ii=1:numel(endog)
            
            subplot(3,2,ii)
            
            d=myirfs.(the_type).(shock).(endog{ii});
            
            out=fanchart(d,ci);
            
            plot_fanchart(out)
            
            title(tex.(endog{ii}))
            
            axis tight
            
        end
        
        [~,h]=sup_label(titel,'t');
        
        set(h,'fontsize',12)
        
    end
    
end