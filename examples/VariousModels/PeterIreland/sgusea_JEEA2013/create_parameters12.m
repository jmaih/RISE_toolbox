function [p,priors]=create_parameters12(ncp_shocks,start_at_mode)
if nargin<2
    start_at_mode=false;
    if nargin<1
        ncp_shocks=false;
    end
end

p=struct();

p.beta = 0.99;
p.mu = 0.34;
p.gam = 2;
p.phi_d = 0.001;
p.delta = 0.025;
p.alpha = 0.36;
p.theta = 1.05;
p.omega = 0.10;

if ncp_shocks
    p.mss_H = 0.9997;
    p.mss_F = 0.9997;
    p.m_hf_ss=1;
    p.zss_H = 1.0045;
else
    p.mss_H = 0.9;
    p.mss_F = 0.9;
    p.zss_H = 1.0042;
end

p.zss_F = p.zss_H;
p.vss_H = 1.0004;
p.vss_F = p.vss_H;
p.z_hf_ss = 1;
p.v_hf_ss = 1;
p.gss_H = 0.20;
p.gss_F = 0.20;
p.rho_m_HF = 0;
p.rho_m_FH = 0;
p.rho_z_HF = 0;
p.rho_z_FH = 0;
p.rho_g_HF = 0;
p.rho_g_FH = 0;
p.rho_v_HF = 0;
p.rho_v_FH = 0;

priors=struct();

if start_at_mode
    priors.theta = {1.5709,0,5};
    
    priors.philtr = {17.113/100,0,3};
    priors.phiktr = {2.971/100,0,3};
    
    if ncp_shocks
        priors.kappamhtr = {0.0077294,0,3};
        priors.kappa_m_F = {8.0501e-06,0,3};
    end
    priors.kappazhtr = {0.002383,0,3};
    priors.kappa_z_F = {1.6159e-05,0,3};
    
    priors.kappavhtr = {0.017416,0,3};
    priors.kappa_v_F = {sqrt(eps),0,3};
    
    priors.rho_m_HH = {0.9671,0,1};
    priors.rho_m_FF = {0.9924,0,1};
    
    priors.rho_z_HH = {0.1752,0,1};
    priors.rho_z_FF = {0.35979,0,1};
    
    priors.rho_v_HH = {0.15793,0,1};
    priors.rho_v_FF = {0.9834,0,1};
    
    priors.rho_g_HH = {0.9695,0,1};
    priors.rho_g_FF = {0.9414,0,1};
    
    priors.sig_m_H = {0.017084,0,2};
    priors.sig_m_F = {0.0063409,0,2};
    
    priors.sig_z_H = {0.01176,0,2};
    priors.sig_z_F = {0.0085657,0,2};
    
    priors.sig_v_H = {0.01436,0,2};
    priors.sig_v_F = {0.00060192,0,2};
    
    priors.sig_g_H = {0.018465,0,2};
    priors.sig_g_F = {0.0098922,0,2};
else
    thetatr = 1.5709;
    priors.theta = {abs(thetatr),0,5};
    
    priors.philtr = {17.1310/100,0,3};
    priors.phiktr = {2.7056/100,0,3};
    
    if ncp_shocks
        kappamftr = 0.0010;
        priors.kappamhtr = {0.0077,0,3};
        priors.kappa_m_F = {abs(kappamftr),0,3};
    end

    priors.kappazhtr = {0.0018,0,3};
    kappazftr = 0.0010;
    priors.kappa_z_F = {abs(kappazftr),0,3};
    
    priors.kappavhtr = {0.0089,0,3};
    kappavftr = 0.0010;
    priors.kappa_v_F = {abs(kappavftr),0,3};
    
    rhomhhtr = sqrt(0.9671/(1-0.9671));
    priors.rho_m_HH = {rhomhhtr^2/(1+rhomhhtr^2),0,1};
    rhomfftr = sqrt(0.9924/(1-0.9924));
    priors.rho_m_FF = {rhomfftr^2/(1+rhomfftr^2),0,1};
    
    rhozhhtr = sqrt(0.1752/(1-0.1752));
    priors.rho_z_HH = {rhozhhtr^2/(1+rhozhhtr^2),0,1};
    rhozfftr = sqrt(0.3598/(1-0.3598));
    priors.rho_z_FF = {rhozfftr^2/(1+rhozfftr^2),0,1};
    
    rhovhhtr = sqrt(0.1581/(1-0.1581));
    priors.rho_v_HH = {rhovhhtr^2/(1+rhovhhtr^2),0,1};
    rhovfftr = sqrt(0.9834/(1-0.9834));
    priors.rho_v_FF = {rhovfftr^2/(1+rhovfftr^2),0,1};
    
    rhoghhtr = sqrt(0.9695/(1-0.9695));
    priors.rho_g_HH = {rhoghhtr^2/(1+rhoghhtr^2),0,1};
    rhogfftr = sqrt(0.9414/(1-0.9414));
    priors.rho_g_FF = {rhogfftr^2/(1+rhogfftr^2),0,1};
    
    sigmhtr = 0.01;
    priors.sig_m_H = {abs(sigmhtr),0,2};
    sigmftr = 0.01;
    priors.sig_m_F = {abs(sigmftr),0,2};
    
    sigzhtr = 0.01;
    priors.sig_z_H = {abs(sigzhtr),0,2};
    sigzftr = 0.01;
    priors.sig_z_F = {abs(sigzftr),0,2};
    
    sigvhtr = 0.01;
    priors.sig_v_H = {abs(sigvhtr),0,2};
    sigvftr = 0.01;
    priors.sig_v_F = {abs(sigvftr),0,2};
    
    sigghtr = 0.01;
    priors.sig_g_H = {abs(sigghtr),0,2};
    siggftr = 0.01;
    priors.sig_g_F = {abs(siggftr),0,2};
end

% add the initial conditions of the priors to the parameters
%------------------------------------------------------------
fields=fieldnames(priors);
for ip=1:numel(fields)
    name=fields{ip};
    p.(name)=priors.(name){1};
end

end