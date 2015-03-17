function [p,priors]=create_parameters(start_at_mode)
if nargin==0
    start_at_mode=false;
end

p=struct();

betatr = sqrt(0.99/0.01);
p.beta = betatr^2/(1+betatr^2);

deltactr = sqrt(0.025/0.975);
p.delta_C = deltactr^2/(1+deltactr^2);

deltaitr = sqrt(0.025/0.975);
p.delta_I = deltaitr^2/(1+deltaitr^2);

priors=struct();
if start_at_mode
    gammatr =   0.622952139588236 ;
    priors.gam = {gammatr^2/(1+gammatr^2),0,1};
    
    thetactr=   0.737470506192051 ;
    priors.theta_C = {thetactr^2/(1+thetactr^2),0,1};
    
    priors.phik_C_tr={0.158646106746036,0,5};
    
    priors.phih_C_tr={0.040927430502772,0,5};
    
    priors.agtr={-0.000046680267058,-2,2};
    
    priors.zg_C_tr={ 0.004778198273099,-2,2};
    
    priors.zg_I_tr={-0.006357571598737,-2,2};
    
    rhoaltr = 0.729388815805754;
    priors.rhola = {rhoaltr^2/(1+rhoaltr^2),0,1};
    
    rhoagtr = -0.000054232696990;
    priors.rhoga = {rhoagtr^2/(1+rhoagtr^2),0,1};
    
    rhocltr =0.909996594639122;
    priors.rhol_C = {rhocltr^2/(1+rhocltr^2),0,1};
    
    rhocgtr = -0.000000834435676*100;
    priors.rhog_C = {rhocgtr^2/(1+rhocgtr^2),0,1};
    
    rhoiltr = 2.118336630027291;
    priors.rhol_I = {rhoiltr^2/(1+rhoiltr^2),0,1};
    
    rhoigtr = 0.647761459082147;
    priors.rhog_I = {rhoigtr^2/(1+rhoigtr^2),0,1};
    
    sigaltr = 0.033038068012394;
    priors.sigla = {abs(sigaltr),0,5};
    
    sigagtr = 0.009363357345790;
    priors.sigga = {abs(sigagtr),0,5};
    
    sigcltr = 0.000000001258949;
    priors.sigl_C = {abs(sigcltr),0,5};
    
    sigcgtr = 0.011647119737011;
    priors.sigg_C = {abs(sigcgtr),0,5};
    
    sigiltr = 0.064514765968648;
    priors.sigl_I = {abs(sigiltr),0,5};
    
    sigigtr = 0.000000008801512;
    priors.sigg_I = {abs(sigigtr),0,5};
    
else
    gammatr = sqrt(0.15/0.85);
    priors.gam = {gammatr^2/(1+gammatr^2),0,1};
    
    thetactr = sqrt(0.30/0.70);
    priors.theta_C = {thetactr^2/(1+thetactr^2),0,1};
    
    priors.phik_C_tr={25/100,0,5};
    
    priors.phih_C_tr={25/100,0,5};
    
    priors.agtr={-0.0001,-2,2};
    
    priors.zg_C_tr={0.0050,-2,2};
    
    priors.zg_I_tr={0.0065,-2,2};
    
    rhoaltr = sqrt(0.50/0.50);
    priors.rhola = {rhoaltr^2/(1+rhoaltr^2),0,1};
    
    rhoagtr = sqrt(0.10/0.90);
    priors.rhoga = {rhoagtr^2/(1+rhoagtr^2),0,1};
    
    rhocltr = sqrt(0.50/0.50);
    priors.rhol_C = {rhocltr^2/(1+rhocltr^2),0,1};
    
    rhocgtr = sqrt(0.10/0.90);
    priors.rhog_C = {rhocgtr^2/(1+rhocgtr^2),0,1};
    
    rhoiltr = sqrt(0.50/0.50);
    priors.rhol_I = {rhoiltr^2/(1+rhoiltr^2),0,1};
    
    rhoigtr = sqrt(0.10/0.90);
    priors.rhog_I = {rhoigtr^2/(1+rhoigtr^2),0,1};
    
    sigaltr = 0.05;
    priors.sigla = {abs(sigaltr),0,5};
    
    sigagtr = 0.05;
    priors.sigga = {abs(sigagtr),0,5};
    
    sigcltr = 0.05;
    priors.sigl_C = {abs(sigcltr),0,5};
    
    sigcgtr = 0.05;
    priors.sigg_C = {abs(sigcgtr),0,5};
    
    sigiltr = 0.05;
    priors.sigl_I = {abs(sigiltr),0,5};
    
    sigigtr = 0.05;
    priors.sigg_I = {abs(sigigtr),0,5};
end
% add the initial conditions of the priors to the parameters
%------------------------------------------------------------
fields=fieldnames(priors);
for ip=1:numel(fields)
    name=fields{ip};
    p.(name)=priors.(name){1};
end

end