%% housekeeping
close all,
clc
format long
%%
Params={
    'h'		      , 0.5  , 0.2  , 'beta'   , 0.9;
    'stilde'      , 1.6  , 0.013, 'uniform', 1  ;
    'gam_tau_l'   , 0.4  , 0.2  , 'gamma'  , 0.9;
    'gam_g_l'     , 1.1  , 0.3  , 'gamma'  , 0.9;
    'rhoa'        , 0.8  , 0.1  , 'beta'   , 0.9;
    'rho_g'	      , 0.8  , 0.1  , 'beta'   , 0.9;
    'rho_tau'	  , 0.5  , 0.2  , 'beta'   , 0.9;
    'rho_z'	      , 0.5  , 0.2  , 'beta'   , 0.9;
    'siga'        , 0.005, 0.01 , 'gamma'  , 0.9;
    'sig_g'	      , 0.02 , 0.015, 'gamma'  , 0.9;
    'sig_tau'	  , 0.005, 0.01 , 'gamma'  , 0.9;
    'sig_z'	      , 0.5  , 0.1  , 'gamma'  , 0.9;
    };
npar=size(Params,1);
bounds=nan(npar,2);
for ii=1:npar
    distr=Params{ii,4};
    m=Params{ii,2};
    s=Params{ii,3};
    prob=Params{ii,5};
    bounds(ii,:)=distributions.find_bounds(distr,m,s,prob);
end
disp(bounds)
%% recover means and standard deviations
MOMS=nan(npar,2);
for ii=1:npar
distr=Params{ii,4};
    prob=Params{ii,5};
    [a,b,moments,fval]=distributions.(distr)(bounds(ii,1),bounds(ii,2),prob);
    MOMS(ii,1)=moments.mean;
    MOMS(ii,2)=moments.sd;
end
%% display and compare
disp('================ Means and Standard Deviations ================ ')
disp(num2cell(MOMS))
disp('================ Discrepancies ================')
disp(abs(cell2mat(Params(:,2:3))-MOMS))
%% housekeeping
format short