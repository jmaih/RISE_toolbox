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

ci=[30,50,68,90];

%% choose identification
quick=true;

Rfunc_star=Rfunc1;

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
        params.(model_types{imod}),Rfunc,shock_names);
    
end

%% plot decompositions

for imod=1:numel(model_types)
    
    the_type=model_types{imod};    

    for iv=1:numel(endog)
        
        vname=tex.(endog{iv});
        
        titel=['Model (',model_types{imod},'): Historical Decomposition of ',vname];
        
        figure('name',titel);
        
        d=pages2struct(hd.(model_types{imod}).(endog{iv}));
        
        contributors=fieldnames(d);
        
        for ii=1:numel(contributors)
            
            subplot(5,2,ii)
            
            out=fanchart(d.(contributors{ii}),ci);
            
            plot_fanchart(out)
            
            title(contributors{ii})
            
            axis tight
            
        end
        
        [~,h]=sup_label(titel,'t');
        
        set(h,'fontsize',12)
        
    end

end