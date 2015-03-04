@#include "emosp_ssmodel.rs"

parameterization
	% fixed parameters
	%------------------
	g		  ,1 ;
	delta	  ,0.025 ;
	eta		  ,1.5 ;
	theta	  ,6 ;
	sig_v	  ,0.01 ;
	% estimated parameters
	%-----------------------
	beta	   ,0.9980 , 0.9, 1;
	gam		   ,0.0736 , 0.0005, 2 ;
	alpha	   ,0.2022 , 0.0005, 1 ;
	phi_p_trans,54.0745/100, 0, 2;
	phi_k_trans,12.4368/100, 0, 2;
	mu_ss_trans,1.0110-1, 0.0005, 2 ;
	omega_r	   ,3.0296 , 0.5, 5;
	omega_mu   ,0.9840 , -3,  3;
	omega_y	   ,-0.0239 , -3,  3 ;
	omega_pai  ,2.0070 , -3,  3 ;
	e_ss	   ,2.7599 , 0.0005, 5  ;
	z_ss_trans ,7034.6/10000, 0, 2;
	rho_a	   ,0.9903 , 0, 1 ;
	rho_e	   ,0.9497 , 0, 1 ;
	rho_x	   ,0.6975 , 0, 1 ;
	rho_z	   ,0.9787 , 0, 1 ;
	rho_v	   ,0.4400 , 0, 1 ;
	sig_a	   ,0.0064 , 0.0005, 2 ;
	sig_e	   ,0.0115 , 0.0005, 2  ;
	sig_x	   ,0.0224 , 0.0005, 2  ;
	sig_z	   ,0.0153 , 0.0005, 2  ;
	