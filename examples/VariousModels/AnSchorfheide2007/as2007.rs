% Sungbae An and Frank Schorfheide (2007): "Bayesian Analysis of DSGE Models"
% Econometric Reviews, 26(2–4):113–172, 2007

endogenous AC "Adjustment costs" B "Bonds" C "Consumption" D "profits" GBAR "Government spending"
G H "Hours worked" LAMBDA M N PAI "Inflation" Q10 R "Interest rate" RSTAR SC T W XI YSTAR Y Z
YGR INFL INT 


exogenous EPS_Z "TFP" EPS_R "Monetary policy" EPS_G "Gov. spending"

parameters upsilon rho_z r_a gam_q kappa tau chi_m chi_h pai_a paistar psi1 psi2 ginv rho_g rho_r
sig_r100 sig_g100 sig_z100

observables YGR INFL INT

model
	% auxiliary parameters
	%----------------------
	# gam = 1 + gam_q/100;
	# beta = 1/(1+r_a/400);
	# pai = 1+pai_a/400;
	# sig_r = sig_r100/100;
	# sig_g = sig_g100/100;
	# sig_z = sig_z100/100;
	# phi = (1-upsilon)*tau/(upsilon*kappa*pai^2); 
	# g = 1/ginv;
	# r = 1+r_a/400; 

	% main model equations
	%----------------------
	Y = N;
	
	AC = phi/2*(PAI-steady_state(PAI))^2;
	
	log(Z) = rho_z*log(Z{-1})+sig_z*EPS_Z;
	
	D = Y - W*N - AC;

	1-phi*PAI*(PAI-steady_state(PAI))-1/upsilon*(1-W-phi/2*(PAI-steady_state(PAI))^2)+
		beta*Q10*phi*(PAI{+1}-steady_state(PAI))*PAI{+1}*gam*Z{+1}*Y{+1}/Y = 0;

	C + B + M - M{-1}/(PAI*gam*Z) + T = W*H + R{-1}*B{-1}/(PAI*gam*Z) + D + SC;

	C^(-tau) - LAMBDA = 0;

	chi_m/M - LAMBDA + beta*LAMBDA{+1}/(PAI{+1}*gam*Z{+1}) = 0;	 

	-chi_h + LAMBDA*W = 0;

	-LAMBDA + beta*R*LAMBDA{+1}/(PAI{+1}*gam*Z{+1}) = 0;

	R = RSTAR^(1-rho_r)*R{-1}^rho_r*exp(sig_r*EPS_R);

	RSTAR = r*paistar*(PAI/paistar)^psi1*(
		@#if gap_rule
			Y/YSTAR 
		@#else
			Z*Y/Y{-1} 
		@#end
		)^psi2;	 

	GBAR + R{-1}*B{-1}/(PAI*gam*Z) = T + B + M - M{-1}/(PAI*gam*Z);

	GBAR = XI*Y;

	G = 1/(1-XI);

	log(G) = (1-rho_g)*log(g)+rho_g*log(G{-1})+sig_g*EPS_G;

	Y = C + GBAR + AC;

	H = N;

	Q10 = (C{+1}/C)^(-tau)*1/(gam*Z{+1});

	YSTAR = G*((1-upsilon)/chi_h)^(1/tau);

	B = 0;

	% Measurement equations x=log(x/x_ss)
	%--------------------------------------

	YGR = gam_q +100*(log(Y/Y{-1})+log(Z/steady_state(Z)));

	INFL = pai_a+400*log(PAI/steady_state(PAI));

	INT = pai_a + r_a + 4*gam_q+400*log(R/steady_state(R));
	

steady_state_model

	paistar=pai;

	G = g;

	Z = 1;

	AC = 0;

	W = 1-upsilon;

	B = 0;
	
	YSTAR = G*((1-upsilon)/chi_h)^(1/tau);

	LAMBDA = chi_h/W;

	C = LAMBDA^(-1/tau);

	XI = 1- 1/G;

	Y = (C+AC)/(1-XI);

	GBAR = XI*Y;

	D = (1-W)*Y - AC;

	N = Y;

	H = N;

	@#if gap_rule
		RSTAR = r*paistar*(Y/YSTAR)^psi2; 
	@#else
		RSTAR = r*paistar*Z^psi2; 
	@#end

	R = RSTAR;

	PAI = beta*R/(gam*Z);

	M = chi_m/((1-beta/(PAI*gam*Z))*LAMBDA);

	T = GBAR-(1-R/(PAI*gam*Z))*B-(1-1/(PAI*gam*Z))*M;

	SC = C + GBAR-W*H-D;

	YGR = gam_q;

	INFL = pai_a;

	INT = pai_a + r_a + 4*gam_q;

parameterization
	chi_h, 1;
	chi_m, 1;
	tau,   2.0, 2.0, 0.5, gamma_pdf;
	kappa, 0.2, 0.2, 0.2, gamma_pdf;
	psi1, 1.5, 1.5, 0.25, gamma_pdf;
	psi2, 0.5, 0.5, 0.25, gamma_pdf;
	rho_r, 0.5, 0.5, 0.2, beta_pdf;
	rho_g, 0.8, 0.8, 0.1, beta_pdf;
	rho_z, 0.66, 0.66, 0.15, beta_pdf;
	r_a,   0.5, 0.5, 0.5, gamma_pdf;
	pai_a, 7.0, 7.0, 2.0, gamma_pdf;
	gam_q, 0.4, 0.4, 0.2, normal_pdf;
	sig_r100, 0.4,0.4, 4, inv_gamma_pdf;
	sig_g100, 1.0, 1.0,4, inv_gamma_pdf;
	sig_z100, 0.5,0.5, 4, inv_gamma_pdf;
	upsilon, 0.1, 0.1, 0.05, beta_pdf;
	ginv, 0.85, 0.85, 0.1, beta_pdf;
