function [p,priors]=create_parameters(version,linear,start_from_mode)
if nargin<3
    start_from_mode=[];
    if nargin<2
        linear=[];
        if nargin<1
            version=[];
        end
    end
end
if isempty(linear)
    error('model type should be specified through linear set to true or false')
end
if isempty(start_from_mode)
    start_from_mode=false;
end
if isempty(version)
    version='a';
end
switch version
    case 'a' % model with an endogenous inflation target.
    case 'b' % model with an exogenous inflation target.
    case 'd' % model with backward-looking price setting.
    otherwise
end

p=struct();
if ~linear
    % own parameterization: those parameters do not exist
    %----------------------------------------------------
    p.thetass=6;
    p.ess=1;
end

% some parameters (beta and zss) are computed from the data
%--------------------------------------------
[~,p]=create_data(p);

% fixed paramters
%-----------------
p.psi = 0.10;
p.delta_e=0;
p.delta_z=0;
p.alpha=1;
rhoxtr = 0/10;
p.rho_x = 10*abs(rhoxtr);
deltaatr = 0;
p.delta_a = abs(deltaatr/10);

% estimated parameters
%-----------------------
priors=struct();
if start_from_mode
    priors.gam ={0.2553,0,1};
    if ~strcmp(version,'d')
        priors.alpha ={0.00001,0,1};
    end
    priors.rho_pai = {0.9069,0,5};
    priors.rho_gy = {0.2347,0,5};
    priors.rho_a = {0.9105,0,1};
    priors.rho_e = {0.0060,0,1};
    priors.rho_v = {0.0546,0,1};
    
    priors.sig_a = {0.0281,0,1};
    priors.sig_e = {0.0007,0,1};
    priors.sig_z = {0.0134,0,1};
    priors.sig_v = {0.0027,0,1};
    priors.sig_pai = {1e-8,0,1};
    if ~strcmp(version,'b')
        priors.delta_e = {0.0010,0,1};
        priors.delta_z = {0.0002,0,1};
    end
else
    gammatr = sqrt(0.50/0.50)/100;
    priors.gam ={(100*gammatr)^2/(1+(100*gammatr)^2),0,1};
    if ~strcmp(version,'d')
        alphatr = sqrt(0.50/0.50)/100;
        priors.alpha ={(100*alphatr)^2/(1+(100*alphatr)^2),0,1};
    end
    rhoptr = 0.25/10;
    priors.rho_pai = {10*abs(rhoptr),0,5};
    rhogytr = 0.05/10;
    priors.rho_gy = {10*abs(rhogytr),0,5};
    
    rhoatr = sqrt(0.75/0.25)/100;
    priors.rho_a = {(100*rhoatr)^2/(1+(100*rhoatr)^2),0,1};
    rhoetr = sqrt(0.50/0.50)/100;
    priors.rho_e = {(100*rhoetr)^2/(1+(100*rhoetr)^2),0,1};
    rhovtr = sqrt(0.25/0.75)/100;
    priors.rho_v = {(100*rhovtr)^2/(1+(100*rhovtr)^2),0,1};
    
    sigatr = 0.01;
    priors.sig_a = {abs(sigatr),0,1};
    sigetr = 0.001;
    priors.sig_e = {abs(sigetr),0,1};
    sigztr = 0.01;
    priors.sig_z = {abs(sigztr),0,1};
    sigvtr = 0.0025;
    priors.sig_v = {abs(sigvtr),0,1};
    sigptr = 0.0005;
    priors.sig_pai = {abs(sigptr),0,1};
    
    if ~strcmp(version,'b')
        deltaetr = 0.01;
        priors.delta_e = {abs(deltaetr/10),0,1};
        deltaztr = 0.001;
        priors.delta_z = {abs(deltaztr/10),0,1};
    end
end

% add the initial conditions of the priors to the parameters
%------------------------------------------------------------
fields=fieldnames(priors);
for ip=1:numel(fields)
    name=fields{ip};
    p.(name)=priors.(name){1};
end

end