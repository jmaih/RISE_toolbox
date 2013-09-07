%% housekeeping
close all
clear
clc
%% add the necessary paths
rise_startup()
%% read the model 
% choose whether to save the generate functions to disk or to keep them in
% memory. I still don't know which one is the quickest...
write_functions_to_disk=false;

test=rise('Canonical_osr');

if write_functions_to_disk
    test=test.set_options('write_functions_to_disk',true);
end

%% estimate with fmincon
tic,test=test.estimate;todisk=toc;

%% speed competition
% test=rise('Canonical_osr');
% test0=rise('Canonical_osr','write_functions_to_disk',true);
% tic,test=test.estimate;todisk=toc;
% tic,test0=test0.estimate;nottodisk=toc;

