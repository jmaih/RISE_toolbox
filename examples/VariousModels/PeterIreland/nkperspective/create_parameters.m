function [p,priors]=create_parameters(start_from_mode)
if nargin==0
    start_from_mode=false;
end

p=struct();

% own parameterization
%----------------------
p.thetass=6;
p.ess=1;
% add zss, paiss and beta
%-------------------------
p.rho_r=1;
p.rho_x=0;
p.psi=0.1;
[~,p]=create_data(p);
% estimated parameters
%----------------------
priors=struct();
if start_from_mode
    priors.gam={0.3904,0,1};
    priors.alpha={0.00001,0,1};
    priors.rho_pai={0.4153,0,3};
    priors.rho_g={0.1270,0,5};
    priors.rho_a={0.9797,0,1};
    priors.rho_e={0.00001,0,1};
    priors.sig_a={0.0868,0,3};
    priors.sig_e={0.0017,0,3};
    priors.sig_z={0.0095,0,3};
    priors.sig_r={0.0014,0,3};
else
    gammatr = sqrt(0.50/0.50);
    priors.gam={gammatr^2/(1+gammatr^2),0,1};
    alphatr = sqrt(0.50/0.50);
    priors.alpha={alphatr^2/(1+alphatr^2),0,1};
    priors.rho_pai={0.25,0,3};
    priors.rho_g={0.05,0,5};
    rhoatr = sqrt(0.75/0.25);
    priors.rho_a={rhoatr^2/(1+rhoatr^2),0,1};
    rhoetr = sqrt(0.25/0.75);
    priors.rho_e={rhoetr^2/(1+rhoetr^2),0,1};
    priors.sig_a={0.01,0,3};
    priors.sig_e={0.001,0,3};
    priors.sig_z={0.01,0,3};
    priors.sig_r={0.0025,0,3};
end

% add the initial conditions of the priors to the parameters
%------------------------------------------------------------
fields=fieldnames(priors);
for ip=1:numel(fields)
    name=fields{ip};
    p.(name)=priors.(name){1};
end

end