%% housekeeping
close all
clear
clc
%% load the data and the estimated models
load tut01_data

load tut02_estimation

load tut04_identification % we use this just to get the shock names

endog=models.ve.endogenous;

%% compute irfs

model_types=fieldnames(models);

myirfs=struct();

for imod=1:numel(model_types)
    
    myirfs.(model_types{imod})=irf(models.(model_types{imod}),shock_names,40);
    
end

%% IRFs comparison

for ishock=1:numel(shock_names)
    
    shock=shock_names{ishock};
    
    titel=['Impulse responses to a ',shock,' shock'];
    
    figure('name',titel);
    
    for ii=1:numel(endog)
        
        subplot(3,2,ii)
        
        d=[];
        
        for imod=1:numel(model_types)
            
            d=[d,myirfs.(model_types{imod}).(shock).(endog{ii})];
            
        end
        
        plot(d,'linewidth',2)
        
        title(tex.(endog{ii}))
        
        if ii==1
            
            legend(model_types)
            
        end
        
    end
    
    [~,h]=sup_label(titel,'t');
    
    set(h,'fontsize',12)
    
end
