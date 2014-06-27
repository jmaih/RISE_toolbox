%% housekeeping
close all
clear
clc
%% add the necessary paths
disp(upper('this file may need updating... contact junior.maih@gmail.com if it does not work'))

setpaths
%% read the generic model (ACCELERATION IS THE DEFAULT MODE)
tc=rise('Canonical_2_MarkovChains');

%% construct the alternative model in which there is no acceleration
tc(2)=tc(1).set_options('accelerate_solver',false);

%% pick a list of variables of interest 
varlist=char('I','PAI','R','Y','ZGDP','ZI','ZPAI','ZY');

%% solve both models and compare them in terms of solution and computing time

t0=zeros(1,2);
for ii=1:2
    tic, tc(ii)=tc(ii).solve; t0(ii)=toc;
    tc(ii).print_solution(varlist)
end
disp('comparison of accelerated(1) vs non-accelerated(2) solutions')
disp(t0/min(t0))

%% check that the solutions are identical
max(max(max(max(abs(tc(1).T-tc(2).T)))))
max(max(max(max(abs(tc(1).R-tc(2).R)))))
