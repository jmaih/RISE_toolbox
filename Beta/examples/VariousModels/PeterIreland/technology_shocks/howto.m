%% housekeeping
close all
clear
clc
%% add the necessary paths
setpaths

%% linearized version of Peter Ireland's model in his techonology shocks paper
% N.B. The models above are the ones Peter Ireland derives with
% micro-foundations. Below is a model he augments in an ad-hock fashion

data=load('gpr.dat');
data=rise_time_series('1948q2',data,char('ln_g','ln_pai','ln_r')); % <-- data=rise_time_series('1948q1',data,char('ln_g','ln_pai','ln_r'));

%  gt = gpr(:,1);
%  pt = gpr(:,2);
%  rt = gpr(:,3);


mm=rise('tshocksnk_loglinearized','data',data,'data_demean',true);
%  gt = gpr(1:127,1);
%  pt = gpr(1:127,2);
%  rt = gpr(1:127,3);

mm(2)=mm(1).set_options('estim_end_date','1979q4');
%   gt = gpr(128:220,1);
%   pt = gpr(128:220,2);
%   rt = gpr(128:220,3);

mm(3)=mm(1).set_options('estim_start_date','1980q1','estim_end_date','2003q1');
%% estimated models and visualize results
for ii=1:3
    mm(ii)=mm(ii).estimate;
end

close all

for ii=1:3
    mm(ii)=mm(ii).irf(16,[],char('ln_g','ln_pai','ln_r','ln_y'));
    mm(ii).print_estimation_results
end
