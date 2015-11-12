function [p,priors]=create_parameters_nonstationary(sticky_price_model,start_at_mode)
if nargin<2
    start_at_mode=false;
    if nargin<1
        sticky_price_model=true;
    end
end

% fixed parameters
%------------------
pfix=struct();
pfix.g =1 ; % the trend will be removed prior to estimation...
deltatr = sqrt(0.025/0.975);
pfix.delta = deltatr^2/(1+deltatr^2); % 0.025 ;
pfix.eta =1.5 ;
thetatr	 =5 ;
pfix.theta = 1 + abs(thetatr);
pfix.sig_v =0.01 ;
% initially set the sticky price parameter to 0 and update it only if
% necessary
%-------------------------------------------------------------------------
pfix.phi_p_trans=0;

% parameters estimated Maximum Likelihood (or rather uniform priors in
% RISE's language
%----------------------------------------------------------------------
priors=struct();
if start_at_mode
    priors.beta	       ={0.9980 , 0.9, 0.9999};
    priors.gam		   ={0.0736 , 0.0005, 2 };
    priors.alpha	   ={0.2022 , 0.0005, 1 };
    if sticky_price_model
        priors.phi_p_trans ={54.0745/100, 0, 2};
    end
    priors.phi_k_trans ={12.4368/100, 0, 2};
    priors.mu_ss_trans ={1.0110-1, 0.0005, 2 };
    priors.omega_r	   ={3.0296 , 0.5, 5};
    priors.omega_mu    ={0.9840 , -5,  5};
    priors.omega_y	   ={-0.0239 , -5,  5 };
    priors.omega_pai   ={2.0070 , -5,  5 };
    priors.e_ss	       ={2.7599 , 0.0005, 5 };
    priors.z_ss_trans  ={7034.6/10000, 0, 2};
    priors.rho_a	   ={0.9903 , 0, 1 };
    priors.rho_e	   ={0.9497 , 0, 1 };
    priors.rho_x	   ={0.6975 , 0, 1 };
    priors.rho_z	   ={0.9787 , 0, 1 };
    priors.rho_v	   ={0.4400 , 0, 1 };
    priors.sig_a	   ={0.0064 , 0.0005, 2 };
    priors.sig_e	   ={0.0115 , 0.0005, 2  };
    priors.sig_x	   ={0.0224 , 0.0005, 2  };
    priors.sig_z	   ={0.0153 , 0.0005, 2  };
else
    betatr = 0.1094;
    priors.beta	       ={(100*betatr)^2/(1+(100*betatr)^2) , 0.9, 0.9999};
    
    priors.gam	       ={0.0428 , 0.0005, 2 };
    
    alphatr = 0.5116;
    priors.alpha	   ={alphatr^2/(1+alphatr^2) , 0.0005, 1 };
    
    if sticky_price_model
        priors.phi_p_trans={1.5125, 0, 2};
    end
    priors.phi_k_trans={0.3426, 0, 2};
    
    priors.mu_ss_trans={0.0086, 0.0005, 2 };
    
    priors.omega_r	   ={2.4903 , 0.5, 5};
    priors.omega_mu    ={1.1670 , -5,  5};
    priors.omega_y	   ={-0.0384 , -5,  5 };
    priors.omega_pai   ={2.8918 , -5,  5 };
    priors.e_ss	       ={3.7927 , 0.0005, 5 };
    priors.z_ss_trans  ={0.7127, 0, 2};
    priors.rho_a	   ={0.9 , 0, 1 };
    priors.rho_e	   ={0.9 , 0, 1 };
    priors.rho_x	   ={0.9 , 0, 1 };
    priors.rho_z	   ={0.9 , 0, 1 };
    priors.rho_v	   ={0.5 , 0, 1 };
    priors.sig_a	   ={0.0200 , 0.0005, 2 };
    priors.sig_e	   ={0.0085 , 0.0005, 2  };
    priors.sig_x	   ={0.3112 , 0.0005, 2  };
    priors.sig_z	   ={0.0201 , 0.0005, 2  };
end

est_list=fieldnames(priors);
p=pfix;
for iname=1:numel(est_list)
    name=est_list{iname};
    p.(name)=priors.(name){1};
end

end