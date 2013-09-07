endogenous
%----------
	Y,   "Output gap",
	PAI, "Inflation",
	I,   "Fed Funds rate"

exogenous
%----------
	EPAI,  "Phil. curve shock",
	EY,    "IS curve shock",
	EI,    "Taylor rule shock"

parameters
%----------
	alpha_pai1, "$\alpha_{\pi,1}$",
	alpha_pai2, "$\alpha_{\pi,2}$",
	alpha_y,    "$\alpha_{y}$",
	c_pai,      "$c_{\pi}$",
	c_y,        "$c_{y}$",
	beta_y1,    "$\beta_{y,1}$",
	beta_y2,    "$\beta_{y,2}$",
	beta_r,     "$\beta_{r}$",
	sig_pai,    "$\sigma_{\pi}$",
	sig_y,      "$\sigma_{y}$",
	sig_i,      "$\sigma_{i}$"
	rho_i,      "$\rho_{i}$",
	gam_y,      "$\gamma_{y}$",
	gam_pai,    "$\gamma_{\pi}$",
	c_i,        "$c_{i}$"

observables
%----------
	I, Y, PAI

model(linear)
%------------
   # alpha_pai = 1/sig_pai;
   # beta_y    = 1/sig_y;
   # gam_i     = 1/sig_i;
   
   alpha_pai*PAI   = c_pai + alpha_pai1*PAI{-1} + alpha_pai2*PAI{-2} +alpha_y*Y{-1}     + EPAI;
   
   beta_y*Y        = c_y   + beta_y1*Y{-1}      + beta_y2*Y{-2} -beta_r*(I{-1}-PAI{-1}) + EY;
   
   gam_i*I         = c_i   + gam_i*rho_i*I{-1}  +gam_i*(1-rho_i)*(gam_y*Y+gam_pai*PAI)  + EI;

%%% read model
%m=rise('svar_constant_test');
%%%
%calib=struct();
%calib.alpha_pai1= 	0.9;
%calib.alpha_pai2= 	0.01  ;
%calib.alpha_y=    	0.1;
%calib.c_pai=      	0  ;
%calib.c_y=        	0  ;
%calib.beta_y1=    	0.9;
%calib.beta_y2=    	0.01  ;
%calib.beta_r=     	0.1;
%calib.sig_pai=   		0.1;
%calib.sig_y=  		0.1;
%calib.sig_i= 			0.1;
%calib.rho_i= 			0.6;
%calib.gam_y= 			0.5;
%calib.gam_pai= 		1.5;
%calib.c_i= 			0  ;
%%%
%m=set(m,'parameters',calib);
%get(m,'parameters')
%m=solve(m)