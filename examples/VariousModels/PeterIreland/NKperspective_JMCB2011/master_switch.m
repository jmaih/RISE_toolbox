%% housekeeping
clear 
close all
clc
%% RISE the model
m=rise('nkp_switch');

%% calibrated parameters and priors

p=struct();
% own parameterization
%----------------------
p.thetass=6;
p.ess=1;
% add zss, paiss and beta
%-------------------------
p.rho_r=1;
p.psi=0.1;
[~,p]=create_data(p);
% estimated parameters
%----------------------
p.gam=0.3904;
p.alpha=0.00001;
p.rho_a=0.9797;
p.rho_e=0.00001;

p.rho_pai_pol_1=0.4153;
p.rho_pai_pol_2=0.4153/4;
p.rho_x_pol_1=0;
p.rho_x_pol_2=0;
p.rho_g_pol_1=0.1270/2;
p.rho_g_pol_2=0.1270;

p.sig_a_vol_1=0.0868/2;
p.sig_a_vol_2=0.0868;
p.sig_e_vol_1=0.0017/2;
p.sig_e_vol_2=0.0017;
p.sig_z_vol_1=0.0095/2;
p.sig_z_vol_2=0.0095;
p.sig_r_vol_1=0.0014/2;
p.sig_r_vol_2=0.0014;

p.pol_tp_1_2=0.01;  
p.pol_tp_2_1=0.1;  
p.vol_tp_1_2=0.1;  
p.vol_tp_2_1=0.2;

m=set(m,'parameters',p);

%% data

[data]=create_data();

%% push into model

m=set(m,'data',data);

%% solve and filter

m=filter(m);
%% Counterfactuals

hd=historical_decomposition(m);

clc

close all

figure('name','Counterfactual output paths ')

contribvars=hd.GHAT.varnames;

% remove all shocks from the contributing list
base_vars=contribvars-m.exogenous.name;

for ishock=1:nshocks
    
    shock=m.exogenous.name{ishock};
    
    % add the shock of interest to the list
    base_plus_shock=[base_vars,shock];
    
    % select the data for the variables above only
    fkst=sum(hd.GHAT(base_plus_shock),2);
    
    % reconstruct the time series
    fkst=ts(hd.GHAT.start,fkst);
    
    subplot(2,2,ishock)
    
    toplot=[data.GHAT+log(p.zss),fkst+log(p.zss)];
    
    % comment out the line below to plot the growth rates instead of the
    % level
    toplot=cumsum(toplot);
    
    plot(toplot,'linewidth',2)
    
    title(m.exogenous.tex_name{ishock})
    
    if ishock==1
        
        legend({'actual','forecast'})
        
    end
    
end
    
%% 
