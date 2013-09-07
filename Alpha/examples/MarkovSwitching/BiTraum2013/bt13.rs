%------------------------------------------------------%
%  RISE version of paper							   %
% " Estimating Fiscal Limits: The Case of Greece"	   %
% by Huixin Bi and Nora Traum (2013)				   %
%------------------------------------------------------%

endogenous	Y "real output", C "consumption", G "govt. purchases"
A "productivity level", N "labor supply", TAU "tax rate",
B "bonds(one-period)", Q "bond price", Z "transfers",
S "debt-to-gdp ratio", BD "post-default liability"

exogenous EZ "transfer shock" EG "govt. purchases shock" ETAU "tax shock" EA "productivity shock"

parameters a rhoa "$\rho _a$", rho_tau "$\rho _{\tau}$",
tau "$\tau$" gam_tau_l "$\gamma ^{\tau,l}$"
b rho_g "$\rho _g$" g gam_g_l "$\gamma ^{g,l}$" rho_z "$\rho _z$" z
h beta "$\beta$" phi "$\phi$", ptilde, phat, stilde, n, b_y	g_y

% The following parameters will be computed endogenously in the steady state file
parameters shat eta1 "$\eta _1$" eta2 "$\eta _2$"

% this is the new syntax for declaring measurement errors
parameters stderr_Y stderr_G stderr_TAU stderr_B stderr_Q

parameters(bt,2) delta "$\delta$" 

parameters(vol,2) siga "$\sigma _a$", sig_tau "$\sigma _{\tau}$", sig_g sig_z "$\sigma _z$" 
% adding some exogenous transition probabilities for the shocks
parameters vol_tp_1_2 vol_tp_2_1

varobs Y G TAU B Q

% throw in some exogenous as well
 EZ ETAU

log_vars BD S N

%%planner_objective{commitment=a,discount=beta}beta*(1-delta)*(C-h*S)/(C-h*B); 

model
	% auxiliary parameters
	# gam_g = (1-rho_g)*gam_g_l*g/b;
	
	# gam_tau = tau*(1-rho_tau)*gam_tau_l/b;

	% endogenous switching probabilities
	! bt_tp_1_2 = exp(eta1+eta2*S)/(1+exp(eta1+eta2*S));  	% (4) p.5 note we use S instead of S{-1} as in the paper
	
	! bt_tp_2_1 = 1-exp(eta1+eta2*S)/(1+exp(eta1+eta2*S));

	% Model equations
	C + G = Y;	% (1) p.4

	Y = A*N;

	A - a = rhoa*(A{-1}-a) + siga*EA;	% (2) p.5

	TAU*Y + B*Q = (1-delta)*B{-1} + G + Z;	% (3) p.5

	S = B/Y;	% 3rd paragraph p.5

	TAU = (1-rho_tau)*tau + rho_tau*TAU{-1} + sig_tau*ETAU + gam_tau*(BD-b);	% (5) p.6

	BD = (1-delta)*B{-1};	% from (3) p.5

	G = (1-rho_g)*g + rho_g*G{-1} + sig_g*EG -gam_g*(BD-b);	% (6) p.6

	Z - z = rho_z*(Z{-1}-z) + sig_z*EZ;	% (7) p.6

	phi*(C-h*C{-1})/(1-N) = A*(1-TAU);	% (10) p.7

	Q = beta*(1-delta)*(C-h*C{-1})/(C{+1}-h*C);	% (11) p.7

steady_state_model
	BD = 3;
	S =4;
	N =5;

parameterization
	ptilde      , 0.3;
	phat	    , 0.999;
	n   	    , 0.25;
	b_y         , 1.095;
	g_y         , 0.181;
	beta	    , 0.99;
	tau		    , 0.33;
	delta(bt,1) , 0;
	delta(bt,2) , 0.075;%0.05
	a           , 1; % start
	b		    , 1; % start
	g		    , 1; % start
	z		    , 1; % start
	phi		    , 1; % start
	% Calibrated Measurement errors on observables
	stderr_Y           , 1; % start
	stderr_G			, 1; % start
	stderr_TAU			, 1; % start
	stderr_B			, 1; % start
	stderr_Q			, 1; % start
	% Estimated parameters
	h		       , 0.13 ,   0.171758483746479 ,   0.828241516252930 , beta_pdf(0.9);
	stilde         , 1.56 ,   1.577483339501605 ,   1.622516660498396 , uniform_pdf  ;
	gam_tau_l      , 0.22 ,   0.136631839674983 ,   0.775365652793273 , gamma_pdf(0.9);
	gam_g_l        , 1.18 ,   0.657219998133462 ,   1.635423871373575 , gamma_pdf(0.9);
	rhoa           , 0.97 ,   0.614610317336182 ,   0.938897071277634 , beta_pdf(0.9);
	rho_g	       , 0.94 ,   0.614610317336182 ,   0.938897071277634 , beta_pdf(0.9);
	rho_tau	       , 0.94 ,   0.171758483746479 ,   0.828241516252930 , beta_pdf(0.9),0,1;
	rho_z	       , 0.58 ,   0.171758483746479 ,   0.828241516252930 , beta_pdf(0.9);
	vol_tp_1_2	   , 0.15 ,   0.01              ,   0.6                              ;     
	vol_tp_2_1	   , 0.15 ,   0.01              ,   0.6                              ;
	siga(vol,1)    , 0.022,   0.000000084371508 ,   0.024202322748895 , gamma_pdf(0.9);
	sig_g(vol,1)   , 0.041,   0.003038737643150 ,   0.049262471293384 , gamma_pdf(0.9);
	sig_tau(vol,1) , 0.026,   0.000000084371508 ,   0.024202322748895 , gamma_pdf(0.9);
	sig_z(vol,1)   , 0.49 ,   0.347642516835017 ,   0.675048065495412 , gamma_pdf(0.9);
	siga(vol,2)    , 0.022,   0.000000084371508 ,   0.024202322748895 , gamma_pdf(0.9);
	sig_g(vol,2)   , 0.041,   0.003038737643150 ,   0.049262471293384 , gamma_pdf(0.9);
	sig_tau(vol,2) , 0.026,   0.000000084371508 ,   0.024202322748895 , gamma_pdf(0.9);
	sig_z(vol,2)   , 0.49 ,   0.347642516835017 ,   0.675048065495412 , gamma_pdf(0.9);

parameter_restrictions
sig_z(vol,2)>sig_z(vol,1);