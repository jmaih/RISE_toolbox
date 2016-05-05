%-------------------------------------------------------------%
%                   Declarations:
%     names are separated by a comma, a space or both
%     "..." are used to describe the preceding variable
%-------------------------------------------------------------%

%Endogenous variables
endogenous	 Y, "Output gap", PAI, "Inflation", I, "Fed Funds rate"

%Exogenous variables
exogenous EPAI,  "Phil. curve shock", EY, "IS curve shock", EI, "Taylor rule shock"

%parameters
parameters alpha_pai1, "$\alpha_{\pi,1}$", alpha_pai2, "$\alpha_{\pi,2}$", alpha_y, "$\alpha_{y}$", c_pai, "$c_{\pi}$",
c_y, "$c_{y}$", beta_y1, "$\beta_{y,1}$", beta_y2, "$\beta_{y,2}$", beta_r, "$\beta_{r}$",

% observable variables
observables I, Y, PAI


model
	% define xx_ssmdef_iary expression to be inserted into the model equations
	% the strategy is to interpret the inverse of the elements below as
	% standard deviations
    # alpha_pai = 1/sig_pai;
    # beta_y    = 1/sig_y;
    # gam_i     = 1/sig_i;
    % Phillips curve
   alpha_pai*PAI   = c_pai + alpha_pai1*PAI{-1} + alpha_pai2*PAI{-2} +alpha_y*Y{-1} + EPAI;
   
   % IS curve
   beta_y*Y = c_y + beta_y1*Y{-1} + beta_y2*Y{-2} -beta_r*(I{-1}-PAI{-1}) + EY;
   
   % Taylor rule
%   gam_i*I = c_i + gam_i*rho_i*I{-1} - gam_i*(1-rho_i)*(gam_y*Y+gam_pai*PAI) + EI;
   gam_i*I = c_i + gam_i*rho_i*I{-1} +gam_i*(1-rho_i)*(gam_y*Y+gam_pai*PAI) + EI;

steady_state_model
	xx_ssmdef_1=alpha_pai-alpha_pai1-alpha_pai2;
	xx_ssmdef_2=gam_i*(1-rho_i);
	xx_ssmdef_3=beta_r*(gam_pai-1);
	xx_ssmdef_4=c_y-beta_r*c_i/xx_ssmdef_2-xx_ssmdef_3*c_pai/xx_ssmdef_1;
	xx_ssmdef_5=beta_y-beta_y1-beta_y2+beta_r*gam_y+xx_ssmdef_3*alpha_y/xx_ssmdef_1;
	
	Y=xx_ssmdef_4/xx_ssmdef_5;
	PAI=c_pai/xx_ssmdef_1+alpha_y*Y/xx_ssmdef_1;
	I=c_i/xx_ssmdef_2+gam_y*Y+gam_pai*PAI;

% the non-policy parameters never switch, they will be controlled by the const markov chain
parameterization
    alpha_pai1, 	0.9, 	0.05, 	1.5, 	gamma_pdf(0.9);
    alpha_pai2, 	0.05  , 	-1  , 	1  , 	normal_pdf(0.9); 
    alpha_y,    	0.1, 	0.05, 	1.5, 	gamma_pdf(0.9);  
    c_pai,      	0  , 	-1  , 	1  , 	normal_pdf(0.9);  
    c_y,        	0  , 	-1  , 	1  , 	normal_pdf(0.9);
    beta_y1,    	0.9, 	0.1 , 	1.5, 	gamma_pdf(0.9);
    beta_y2,    	0.05  , 	-2  , 	2  , 	normal_pdf(0.9);
    beta_r,     	0.1, 	0.05, 	1  , 	gamma_pdf(0.9); 

	