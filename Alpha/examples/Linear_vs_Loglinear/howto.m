%% housekeeping
% for illustration, we use a model by Peter Ireland
close all
clear
clc
%% add the necessary paths
disp(upper('this file may need updating... contact junior.maih@gmail.com if it does not work'))

setpaths

%% read models with their respective steady state files
pni=rise('tshocksnk','steady_state_file','tshocksnk_steadystate');

pni(2)=rise('tshocksnk_log_linear','steady_state_file','tshocksnk_steadystate_loglinear');

% we loglinearize all the endogenous variables
pni(3)=pni(1).set_options('loglinearized',{pni(1).varendo.name});

%% solve the models, print solutions and compare
for ii=1:3
    pni(ii)=pni(ii).solve;
    pni(ii).print_solution
end


%% look at some impulse responses
pni(1)=pni(1).irf(20,[],char('g','pai','r','y'));

pni(2)=pni(2).irf(20,[],char('ln_g','ln_pai','ln_r','ln_y'));

%% compare the derivatives
% here we show that multiplying the structural matrices of the the
% linearized model by the steady state gives the structural matrices of the
% log-linearized model. Just as the algebra predicts
% let z=exp(x). The f(z)=f(exp(x))=g(x). Taking the Taylor expansion of
% both sides f(z0)+(z-z0)*f'(z0)=g(x0)+(x-x0)*g'(x0). We also have by the
% chain rule that exp(x)*f'(exp(x))=g'(x). Hence
% f(z0)+(z-z0)*f'(z0)=g(x0)+(x-x0)*exp(x0)*f'(exp(x0)) or
% f(z0)+(z-z0)*f'(z0)=g(x0)+(x-x0)*z0*f'(z0). We can express everything in
% terms of z. f(z0)+(z-z0)*f'(z0)=f(z0)+(ln(z)-ln(z0))*z0*f'(z0).
% simplifying, (z-z0)=(ln(z)-ln(z0))*z0 or (z-z0)/z0=(ln(z)-ln(z0))
max(max(abs(bsxfun(@times,pni(1).Aplus,[pni(1).varendo.det_steady_state])-pni(2).Aplus)))
max(max(abs(bsxfun(@times,pni(1).A0,[pni(1).varendo.det_steady_state])-pni(2).A0)))
max(max(abs(bsxfun(@times,pni(1).Aminus,[pni(1).varendo.det_steady_state])-pni(2).Aminus)))
%% linearized version of the model with parameters computed outside the model
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
