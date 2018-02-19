%% housekeeping
close all
clear
clc
%% load the data

load tut01_data

load tut02_estimation

load tut09_bootstrap

load tut04_identification

ci=[30,50,68,90];

%% pick a model

model='ve_lr'; % ve

endog=models.(model).endogenous;

%% compute decompositions

% choose identification scheme
Rfunc=[];

hd=variance_decomposition(models.(model),params.(model),Rfunc,shock_names);

%% plot decompositions
shock_tex=shock_names;

myrange='1:50';

for iv=1:numel(endog)
    
    vname=tex.(endog{iv});
    
    titel=['Variance Decomposition (in %) of ',vname];
    
    figure('name',titel);
    
    d=hd.conditional.(endog{iv});
    
    d.varnames=shock_names;
    
    d=pages2struct(d);
    
    contributors=fieldnames(d); % = shock_names
    
    for ii=1:numel(contributors)
    
        subplot(3,2,ii)
        
        % note we are multiplying by 100, this just by pure convenience
        %--------------------------------------------------------------
        out=fanchart(100*d.(contributors{ii})(myrange),ci);
        
        plot_fanchart(out)
        
        title(contributors{ii})
        
        axis tight
        
    end
    
    [~,h]=sup_label(titel,'t');
    
    set(h,'fontsize',12)
    
end