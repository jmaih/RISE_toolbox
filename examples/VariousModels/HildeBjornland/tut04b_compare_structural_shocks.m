%% housekeeping
close all
clear
clc
%% load the data

load tut02_estimation

load tut04_identification

endog=models.ve.endogenous;

%% choose identification
quick=true;

Rfunc_star=Rfunc1;

params=[];

%% compute shocks

model_types=fieldnames(models);

sshocks=struct();

for imod=1:numel(model_types)
    
    if quick
        
        Rfunc=[];
        
    else
        
        Rfunc=Rfunc_star.(model_types{imod});
        
    end
    
    sshocks.(model_types{imod})=structural_shocks(models.(model_types{imod}),...
        params,Rfunc,shock_names);
    
end

%% plots

titel='Structural shocks';

figure('name',titel);

for ii=1:numel(shock_names)
    
    thisname=shock_names{ii};
    
    subplot(3,2,ii)
    
    d=[];
    
    for imod=1:numel(model_types)
        
        the_type=model_types{imod};
        
        d=[d,sshocks.(the_type).(thisname)];
        
    end
    
    plot(d,'linewidth',2)
    
    title(thisname)
    
    if ii==1
        
        legend(model_types)
        
    end
    
end

[~,h]=sup_label(titel,'t');

set(h,'fontsize',12)


