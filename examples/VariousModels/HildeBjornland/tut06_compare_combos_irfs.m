%% housekeeping
close all
clear
clc
%% load the data

load tut01_data

load tut02_estimation

load tut04_identification

%% IRFs

model_types=fieldnames(models);

endog=models.(model_types{1}).endogenous;

myirfs1=struct();

myirfs2=struct();

myirfs3=struct();

for imod=1:numel(model_types)
    
    % replicating the choleski restrictions
    %--------------------------------------
    myirfs1.(model_types{imod})=irf(models.(model_types{imod}),shock_names,...
        40,[],Rfunc1.(model_types{imod}));
    
    % Choleski (top) + Long run restriksjon på valutakurs
    %-----------------------------------------------------------
    myirfs2.(model_types{imod})=irf(models.(model_types{imod}),shock_names,...
        40,[],Rfunc2.(model_types{imod}));
    
    % Choleski (top) + sign på valutakurs
    %--------------------------------------
    myirfs3.(model_types{imod})=irf(models.(model_types{imod}),shock_names,...
        40,[],Rfunc3.(model_types{imod}));
    
end

%% plots

for imod=1:numel(model_types)
    
    for ishock=1:numel(shock_names)
        
        shock=shock_names{ishock};
        
        titel=['Model (',model_types{imod},'): Impulse responses to a ',shock,' shock'];
        
        figure('name',titel);
        
        for ii=1:numel(endog)
            
            subplot(3,2,ii)
            
            d=[myirfs1.(model_types{imod}).(shock).(endog{ii}),...
                myirfs2.(model_types{imod}).(shock).(endog{ii}),...
                myirfs3.(model_types{imod}).(shock).(endog{ii})];
            
            texname=tex.(endog{ii});
            
            plot(d,'linewidth',2)
            
            title(texname)
            
            if ii==1
                
                legend('choleski','long-run(RER)','sign(RER)')
                
            end
            
        end
        
        [~,h]=sup_label(titel,'t');
        
        set(h,'fontsize',12)
        
    end
    
end
