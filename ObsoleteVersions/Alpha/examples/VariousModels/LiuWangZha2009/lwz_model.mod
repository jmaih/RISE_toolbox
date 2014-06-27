//Liu, Wang, and Zha: Land-price dynamics and macroeconomic fluctuations.
//
//-------------------------------------
//Shocks identified in the model
//-------------------------------------
// Eps_a:         Intertemporal preference (risk premium) shock
// Eps_phi:       Housing demand shock
// Eps_psi:       Labor supply shock
// Eps_z:         Permanent neutral technology shock
// Eps_v:         Transitory neutral technology shock
// Eps_q:         Permanent shock to investment specific technical change
// Eps_mu:        Transitory shock to investment specific technical change
// Eps_xitheta:   Collateral shock


//----------------------------------------------
//--  Declare varaibles, shocks,  and parameter 
//----------------------------------------------

var 
//--- 19 Endogenous state variables
Muh W Ql R Mue Mub N I_m Y Ch Ce Qk Lh Le K B G_Lambda G_Z G_Q
//--- 8 variables appearing in shock pricesses
Lambda_z V Lambda_q Mu_m Lambda_a Phi Psi_m Thetat
//--- 6 Observed data
DLogQl DLogQ  DLogC  DLogI  DLogB  LogL
//--- Additional (auxillary) variables
C DLogGDP
;

varexo Eps_a Eps_phi Eps_psi Eps_z Eps_v Eps_q  Eps_mu  Eps_xitheta;

parameters
//--- 5 estimated parameters
gamma_h gamma_e omega g_trans lambda_qbar_trans
//--- 8 calibrated parameters
alpha_m r ky ik qlLeOY theta_bar qlLhOY n lbar
//--- 8 parameters related to shock processes
rho_z rho_v rho_q rho_mu rho_a rho_phi rho_psi rho_xitheta
// standard deviations
sig_a,      
sig_z,      
sig_v,      
sig_q,      
sig_mu,     
sig_phi,    
sig_psi,    
sig_xitheta
;


//--------------------------------------------------------------------
//-- Model specification
//--------------------------------------------------------------------
model(linear);
	# g = g_trans/100 + 1;
	# lambda_qbar = 0.01*lambda_qbar_trans + 1;
	# lambda_kbar = g*lambda_qbar; //See second sentence of first full paragraph of page 9
	# delta = 1-(1-ik)*lambda_kbar ;  //Rearrange equation 60 
	
	//From Tao's calibration5.m file
	# betaparam = (qlLeOY*(1.0-g*theta_bar/r)+ky*(1.0-g*theta_bar/(lambda_kbar*r)))/(alpha_m+qlLeOY*(1.0-theta_bar)+ky*(1.0-delta-theta_bar)/lambda_kbar); 
	# lambda_abar = g/(betaparam*r)-1.0; //See equation 58 (first part)
	# muboe = (betaparam*lambda_abar)/g; //See equation 58 (second part)
	# phi = qlLeOY*(1.0-betaparam*(1.0+lambda_abar*theta_bar))/(betaparam*alpha_m); // from calibration5.m
	// # qlLeOY = (betaparam*alpha_m*phi)/(1 - betaparam - betaparam*lambda_abar*theta_bar) //Equation 59;
	//Do not need equation 61
	
	# iy = ik*ky; //Equation 62
	# by = theta_bar*g*qlLeOY + (theta_bar/lambda_qbar)*ky; //Equation 63
	
	
	# ceOY = alpha_m-iy-by*(1-betaparam*(1+lambda_abar))/g; //See equation 64
	# chOY = 1-ceOY-iy; //See equation 65
	
	# lhOle = qlLhOY/qlLeOY; //See equation 67
	# le_m = lbar/(1.0+lhOle);
	# lh = lbar-le_m;
	
	# psi_bar = ((1-alpha_m)*g*(1-gamma_h/r)*(1/chOY))/(n*(g-gamma_h)); //Equation 68 
	
	# omegah = (g - betaparam*(1+lambda_abar)*gamma_h)*(g - gamma_h); //See top of page 12
	# omegae = (g - betaparam*gamma_e)*(g - gamma_e); //See top of page 12
	
	//(1)
	omegah*Muh = -(g^2 + (gamma_h^2)*betaparam*(1+lambda_abar))*Ch + g*gamma_h*(Ch(-1) - G_Lambda) -
	betaparam*lambda_abar*gamma_h*(g - gamma_h)*Lambda_a(+1) + 
	betaparam*(1+lambda_abar)*g*gamma_h*(Ch(+1) + G_Lambda(+1)); //Eqn (69)
	
	
	//(2)
	W + Muh = Psi_m ; // Eqn. 70
	
	//(3)
	Ql + Muh = betaparam*(1+lambda_abar)*(Muh(+1) + Ql(+1)) + (1 - betaparam*(1+lambda_abar))*(Phi -
	Lh) + betaparam*lambda_abar*Lambda_a(+1); // Eqn. 71
	
	
	//(4) 
	Muh - R = Muh(+1) + (lambda_abar/(1+lambda_abar))*Lambda_a(+1) - G_Lambda(+1); // Eqn (72)
	
	//(5)
	omegae*Mue = -(g^2 + betaparam*gamma_e^2)*Ce + g*gamma_e*(Ce(-1) - G_Lambda) +
	betaparam*g*gamma_e*(Ce(+1) + G_Lambda(+1)); // Eqn (73)
	
	//(6)
	W = Y - N; // Eqn (74)
	
	//(7)
	Qk = (1+betaparam)*omega*(lambda_kbar^2)*I_m - omega*(lambda_kbar^2)*I_m(-1) +
	omega*(lambda_kbar^2)*(G_Lambda + G_Q) - betaparam*omega*(lambda_kbar^2)*(I_m(+1)
	 + G_Lambda(+1) + G_Q(+1)); // Eqn (75)
	
	//(8)
	Qk + Mue = (muboe*theta_bar/lambda_qbar)*(Mub + Thetat) + 
	(betaparam*(1-delta)/(lambda_kbar))*(Qk(+1) - G_Q(+1) - G_Lambda(+1)) + 
	(1 - (muboe*theta_bar/lambda_qbar))*Mue(+1) + 
	(muboe*theta_bar/lambda_qbar)*(Qk(+1) - G_Q(+1)) + 
	betaparam*alpha_m*(1 - phi)*(1/ky)*(Y(+1) - K); // Eqn (76)
	
	//(9)
	Ql + Mue = muboe*g*theta_bar*(Thetat + Mub) + (1 - muboe*g*theta_bar)*Mue(+1) + 
	muboe*g*theta_bar*(Ql(+1) + G_Lambda(+1)) + betaparam*Ql(+1) + 
	(1 - betaparam - betaparam*lambda_abar*theta_bar)*(Y(+1) - Le); // Eqn (77)
	
	//(10)
	Mue - R = (1/(1+lambda_abar))*(Mue(+1) - G_Lambda(+1) + lambda_abar*Mub); // Eqn (78)
	
	//(11)
	Y = alpha_m*phi*Le(-1) + alpha_m*(1- phi)*K(-1) + (1-alpha_m)*N - 
	(((1-phi)*alpha_m)/(1-((1-phi)*alpha_m)))*(G_Z + G_Q); // Eqn (79)
	
	//(12)
	K = ((1-delta)/(lambda_kbar))*(K(-1) - G_Lambda - G_Q) 
	+ (1 - ((1-delta)/(lambda_kbar)))*I_m; // Eqn (80)
	
	//(13)
	Y = chOY*Ch + ceOY*Ce + iy*I_m; // Eqn (81)
	
	//(14)
	0 = (lh/lbar)*Lh + (le_m/lbar)*Le;  // Eqn (82)
	
	//(15)
	alpha_m*Y = ceOY*Ce + iy*I_m + qlLeOY*(Le - Le(-1)) + ((1/g)*by)*(B(-1) - G_Lambda) 
	- ((1/r)*by)*(B - R);// Eqn (83)
	
	//(16)
	B = Thetat + (g*theta_bar*(qlLeOY/by))*(Ql(+1) + Le + G_Lambda(+1)) + 
	(1 - (g*theta_bar*(qlLeOY/by)))*(Qk(+1) + K - G_Q(+1)); // Eqn (84)
	
	//(17)
	G_Z = Lambda_z + V - V(-1); // Eqn (85)

	//(18)
	G_Q = Lambda_q + Mu_m - Mu_m(-1); // Eqn (86)

	//(19)
	G_Lambda = (1/(1-(1-phi)*alpha_m))*G_Z + (((1-phi)*alpha_m)/(1-(1-phi)*alpha_m))*G_Q; // Eqn (87)
	
	
	//Shock processes
	
	Lambda_z = rho_z*Lambda_z(-1) + sig_z*Eps_z; // Eqn (88)
	
	V = rho_v*V(-1) + sig_v*Eps_v; // Eqn (89)
	
	Lambda_q = rho_q*Lambda_q(-1) + sig_q*Eps_q; // Eqn (90)
	
	Mu_m = rho_mu*Mu_m(-1) + sig_mu*Eps_mu; // Eqn (91)
	
	Lambda_a = rho_a*Lambda_a(-1) + sig_a*Eps_a; // Eqn (93)
	
	Phi = rho_phi*Phi(-1) + sig_phi*Eps_phi; // Eqn (94)
	
	Psi_m = rho_psi*Psi_m(-1) + sig_psi*Eps_psi; // Eqn (95)
	
	Thetat = rho_xitheta*Thetat(-1) + sig_xitheta*Eps_xitheta; //Eqn (96)
	
	C = (chOY/(chOY+ceOY))*Ch + (ceOY/(chOY+ceOY))*Ce;
	
	//--- Observer equations
	DLogQl = (Ql - Ql(-1) + log(g) + G_Lambda); // relative price of liquid land
	DLogQ =  (G_Q + log(lambda_qbar)); // inverse of the relative price of investment
	DLogC = ((chOY/(chOY+ceOY))*(Ch - Ch(-1)) + (ceOY/(chOY+ceOY))*(Ce - Ce(-1)) + log(g) + G_Lambda); //real consumption
	DLogI = (I_m - I_m(-1) + log(g) + G_Lambda); //investment plus durables deflated by consumption deflator
	DLogB = (B - B(-1) + log(g) + G_Lambda); //real debt growth
	LogL = N + log(n); //log hours
	
	//--- Auxillary variables
	DLogGDP = iy*DLogI + (1-iy)*DLogC;
end;


//--------------------------------------------------------------------
//-- Estimation
//--------------------------------------------------------------------
varobs DLogQl DLogQ DLogC DLogI DLogB LogL;

parameterization;
	//--- Targeted steady state values
	alpha_m,            0.30; //(1) the average labor income share is 70% (page 15 of paper)          
	r,                  1.01; //(2) the average real prime loan rate is 4% per annum (page 15 of paper)
	ky,               4.6194; //(3) the capital-output ratio is on average 1.15 at the annual frequency (page 15 of paper)
	ik,             0.2093/4; // (4) the investment-capital ratio is on average 0.209 at the annual frequency  
	qlLeOY,             2.60; //(5) the average land-output ratio is 0.65 at the annual frequency
	theta_bar,          0.75;  //(6) the average nonfarm and nonfiancial businesses' loan-asset ratio is 0.75 at the annual frequency (page 15 of paper)
	qlLhOY,           5.8011; //(7) the average housing output ratio is 1.45 at the annual frequency
	n,                   1/4; //(8) average market hours is 25% of time endowment
	lbar,                  1; //arbitrary
	// PARAM NAME, INITVAL, LB, UB, PRIOR_SHAPE, PRIOR_P1, PRIOR_P2, PRIOR_P3, PRIOR_P4, JSCALE
	// PRIOR_SHAPE: BETA_PDF, GAMMA_PDF, NORMAL_PDF, INV_GAMMA_PDF
	sig_a,             0.1387,	0.0001,   2.0000,    inv_gamma_pdf,  .90;
	sig_z,             0.0036,	0.0001,   2.0000,    inv_gamma_pdf,  .90;
	sig_v,             0.0038,	0.0001,   2.0000,    inv_gamma_pdf,  .90;
	sig_q,             0.0037,	0.0001,   2.0000,    inv_gamma_pdf,  .90;
	sig_mu,            0.0025,	0.0001,   2.0000,    inv_gamma_pdf,  .90;
	sig_phi,           0.0543,	0.0001,   2.0000,    inv_gamma_pdf,  .90;
	sig_psi,           0.0073,	0.0001,   2.0000,    inv_gamma_pdf,  .90;
	sig_xitheta,       0.0126,	0.0001,   2.0000,    inv_gamma_pdf,  .90;
	gamma_h,           0.7000,	0.0250,   0.7760,    beta_pdf,       .90;
	gamma_e,           0.7000,	0.0250,   0.7760,    beta_pdf,       .90;
	omega,             3.0000, 	0.1020,   5.9940,    gamma_pdf,      .90;
	g_trans,           0.3750,	.10000,   1.5000,    gamma_pdf,      .90;
	lambda_qbar_trans, 1.2500,	.10000,   1.5000,    gamma_pdf,      .90;
	rho_a,             0.9000,	0.0250,   0.7760,    beta_pdf,       .90;
	rho_z,             0.9000,	0.0250,   0.7760,    beta_pdf,       .90;
	rho_v,             0.9000,	0.0250,   0.7760,    beta_pdf,       .90;
	rho_q,             0.9000,	0.0250,   0.7760,    beta_pdf,       .90;
	rho_mu,            0.9000,	0.0250,   0.7760,    beta_pdf,       .90;
	rho_phi,           0.9000,	0.0250,   0.7760,    beta_pdf,       .90;
	rho_psi,           0.9000,	0.0250,   0.7760,    beta_pdf,       .90;
	rho_xitheta,       0.9000,	0.0250,   0.7760,    beta_pdf,       .90;
end;
