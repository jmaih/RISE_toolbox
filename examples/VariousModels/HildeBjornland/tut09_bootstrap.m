%% housekeeping
close all
clear
clc
%% load the data

load tut02_estimation

model_types=fieldnames(models);

%% bootstrap

n=1000;

params=struct();

for ii=1:numel(model_types)
    
    params.(model_types{ii})=bootstrap(models.(model_types{ii}),n);
    
end

%% save sampled parameters for later use

save('tut09_bootstrap','params')
									
