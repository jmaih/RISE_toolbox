%% housekeeping
close all
clear
clc
%% load the data

load tut01_data

load tut02_estimation

load tut04_identification

endog=models.ve.endogenous;

ci=[30,50,68,90];

%% choose identification
quick=true;

Rfunc_star=Rfunc1;

params=[];

%% compute decompositions

model_types=fieldnames(models);

hd=struct();

for imod=1:numel(model_types)
    
    if quick
        
        Rfunc=[];
        
    else
        
        Rfunc=Rfunc_star.(model_types{imod});
        
    end
    
    hd.(model_types{imod})=historical_decomposition(models.(model_types{imod}),...
        params,Rfunc,shock_names);
    
end

%% plot decompositions
shock_only=true;

for imod=1:numel(model_types)
    
    the_type=model_types{imod};
    
    if strcmp(the_type,'ve')
        
        inBetween='out';
        
    else
        
        inBetween='';
        
    end
    
    titel=['Model with',inBetween,' block exogeneity: Historical Decomposition'];
    
    figure('name',titel);
    
    for iv=1:numel(endog)
        
        d=hd.(model_types{imod}).(endog{iv});
        
        subplot(3,2,iv)
        
        if shock_only
            
            shock_tex=shock_names;
            
            d=d(shock_names);
            
        else
            
            shock_tex=d.varnames;
            
        end
        
        plot_decomp(d)
        
        if iv==1
            
            legend(shock_tex,'location','SE',...
                'Orientation','horizontal')
            
        end
        
        title(tex.(endog{iv}))
        
    end
    
    [~,h]=sup_label(titel,'t');
    
    set(h,'fontsize',12)
    
end