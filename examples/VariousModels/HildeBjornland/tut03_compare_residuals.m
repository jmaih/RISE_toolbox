%% housekeeping
close all
clear
clc
%% load the data

load tut02_estimation

params=[];

%% compute residuals

re=struct();

model_types=fieldnames(models);

for imod=1:numel(model_types)
    
    re.(model_types{imod})=residuals(models.(model_types{imod}),params);
    
end

%% plot residuals

res_names=fieldnames(re.(model_types{1}));

titel='Residuals';

figure('name',titel);

for ii=1:numel(res_names)
    
    thisname=res_names{ii};
    
    subplot(3,2,ii)
    
    d=[];
    
    for imod=1:numel(model_types)
        
        d=[d,re.(model_types{imod}).(thisname)];
        
    end
    
    plot(d,'linewidth',2)
    
    title(thisname)
    
end
    
    [~,h]=sup_label(titel,'t');
    
    set(h,'fontsize',12)

