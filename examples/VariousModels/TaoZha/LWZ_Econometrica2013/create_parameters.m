function [p,priors]=create_parameters(start_from_mode)
p=struct();
%--- Targeted steady state values
p.alpha = 0.30; %(1) the average labor income share is 70% (page 15 of paper)
p.r = 1.01; %(2) the average real prime loan rate is 4% per annum (page 15 of paper)
p.ky = 4.6194; %(3) the capital-output ratio is on average 1.15 at the annual frequency (page 15 of paper)
p.ik = 0.2093/4; % (4) the investment-capital ratio is on average 0.209 at the annual frequency
p.qlLeOY = 2.60; %(5) the average land-output ratio is 0.65 at the annual frequency
p.theta_bar = 0.75;  %(6) the average nonfarm and nonfiancial businesses' loan-asset ratio is 0.75 at the annual frequency (page 15 of paper)
p.qlLhOY = 5.8011; %(7) the average housing output ratio is 1.45 at the annual frequency
p.n = 1/4; %(8) average market hours is 25% of time endowment
p.lbar = 1; %arbitrary

priors=struct();
priors.gamma_h          ={ 0.7  , 0.0256, 0.7761 , 'beta_pdf(0.9)'     , .00000001     , .99999999999};
priors.gamma_e          ={ 0.7  , 0.0256, 0.7761 , 'beta_pdf(0.9)'     , .00000001     , .99999999999};
priors.omega            ={ 3.0  , 0.1020, 5.994  , 'gamma_pdf(0.9)'    , .00001        , 10          };
priors.g_trans          ={ 0.375, 0.1   , 1.5    , 'gamma_pdf(0.9)'    , .00001        , 10          };
priors.lambda_qbar_trans={ 1.25 , 0.1   , 1.5    , 'gamma_pdf(0.9)'    , .00001        , 10          };
priors.rho_a            ={ 0.9  , 0.0256, 0.7761 , 'beta_pdf(0.9)'     , .00000001     , .99999999999};
priors.rho_z            ={ 0.9  , 0.0256, 0.7761 , 'beta_pdf(0.9)'     , .00000001     , .99999999999};
priors.rho_v            ={ 0.9  , 0.0256, 0.7761 , 'beta_pdf(0.9)'     , .00000001     , .99999999999};
priors.rho_q            ={ 0.9  , 0.0256, 0.7761 , 'beta_pdf(0.9)'     , .00000001     , .99999999999};
priors.rho_mu           ={ 0.9  , 0.0256, 0.7761 , 'beta_pdf(0.9)'     , .00000001     , .99999999999};
priors.rho_phi          ={ 0.9  , 0.0256, 0.7761 , 'beta_pdf(0.9)'     , .00000001     , .99999999999};
priors.rho_psi          ={ 0.9  , 0.0256, 0.7761 , 'beta_pdf(0.9)'     , .00000001     , .99999999999};
priors.rho_xitheta      ={ 0.9  , 0.0256, 0.7761 , 'beta_pdf(0.9)'     , .00000001     , .99999999999};
priors.sig_Eps_a        ={ 0.01 , 0.0001, 2      , 'inv_gamma_pdf(0.9)', 0.000000000001, 100         };
priors.sig_Eps_z        ={ 0.01 , 0.0001, 2      , 'inv_gamma_pdf(0.9)', 0.000000000001, 100         };
priors.sig_Eps_v        ={ 0.01 , 0.0001, 2      , 'inv_gamma_pdf(0.9)', 0.000000000001, 100         };
priors.sig_Eps_q        ={ 0.01 , 0.0001, 2      , 'inv_gamma_pdf(0.9)', 0.000000000001, 100         };
priors.sig_Eps_mu       ={ 0.01 , 0.0001, 2      , 'inv_gamma_pdf(0.9)', 0.000000000001, 100         };
priors.sig_Eps_phi      ={ 0.01 , 0.0001, 2      , 'inv_gamma_pdf(0.9)', 0.000000000001, 100         };
priors.sig_Eps_psi      ={ 0.01 , 0.0001, 2      , 'inv_gamma_pdf(0.9)', 0.000000000001, 100         };
priors.sig_Eps_xitheta  ={ 0.01 , 0.0001, 2      , 'inv_gamma_pdf(0.9)', 0.000000000001, 100         };

if start_from_mode
    priors.gamma_e{1}= 0.658435028;
    priors.gamma_h{1}= 0.497555761;
    priors.omega{1}= 0.175347074;
    priors.g_trans{1}= 0.422109209;
    priors.lambda_qbar_trans{1}= 1.21261837;
    priors.rho_a{1}=  0.905471687;
    priors.rho_v{1}= 0.00947706131;	%rho_{nu_z}
    priors.rho_mu{1}= 0.294889827;	   %rho_{nu_q}
    priors.rho_psi{1}=  0.982918278;
    priors.rho_q{1}= 0.561965864;
    priors.rho_xitheta{1}= 0.980427788;
    priors.rho_phi{1}=  0.999758433;
    priors.rho_z{1}= 0.426298781;	 %rho_z
    priors.sig_Eps_z{1}=     0.00418908791;
    priors.sig_Eps_v{1}=   0.00365910563;
    priors.sig_Eps_q{1}=     0.00415238812;
    priors.sig_Eps_mu{1}=   0.00288990925;
    priors.sig_Eps_a{1}=     0.101281626;
    priors.sig_Eps_phi{1}=   0.0462334747;
    priors.sig_Eps_psi{1}=    0.00731583619;
    priors.sig_Eps_xitheta{1}= 0.0111777605;
end

% add the initial conditions of the priors to the parameters
%------------------------------------------------------------
fields=fieldnames(priors);
for ip=1:numel(fields)
    name=fields{ip};
    p.(name)=priors.(name){1};
end

end