%---------------------------------------------------------------------------------
% Liu, Wang, and Zha(2013): "Land-price dynamics and macroeconomic fluctuations."
%  Econometrica, Vol. 81, No. 3, May, pages 1147-1184
% Fval obtained by the minimization routine: -2430.913404
%---------------------------------------------------------------------------------

endogenous 
%--- 19 Endogenous state variables
Muh "Lagrange mult. on flow of funds(household)"
W "Wage rate"
Ql "Relative price of land"
R "Gross real loan rate"
Mue "Lagrange mult. on flow of funds(entrepreneur)"
Mub "Lagrange mult. on borrowing constr."
N "Labor hours"
I "Investment"
Y "Output"
Ch "Consumption (household)"
Ce "Consumption(entrepreneur)"
Qk
Lh "Land holdings (household)"
Le "Land holdings(entrepreneur)"
K "Capital"
B "Debt"
G_Lambda "Growth rate of composite stochastic trend" % G_gam in text
G_Z
G_Q "Growth rate of stochastic process Q"
%--- 8 variables appearing in shock processes
Lambda_z "Perm. neutral techn. process"
V "Transitory neutral tech. process"
Lambda_q "Perm. investment specific techn. change process"
Mu "Transitory investment specific techn. change process"
Lambda_a "Intertemporal preference (risk premium) process"
Phi "Housing demand process"
Psi "Labor supply process"
Thetat "Collateral process"
%--- 6 Observed data
DLogQl "relative price of liquid land"
DLogQ "inverse of the relative price of investment"
DLogC "real consumption"
DLogI "investment + durables deflated by consumption deflator"
DLogB "real debt growth"
LogL "log hours"
%--- Additional (auxillary) variables
C "Total consumption" DLogGDP

observables DLogQl DLogQ DLogC DLogI DLogB LogL

exogenous
Eps_a "Intertemporal preference (risk premium) shock"
Eps_phi "Housing demand shock"
Eps_psi "Labor supply shock"
Eps_z "Permanent neutral technology shock"
Eps_v "Transitory neutral technology shock"
Eps_q "Permanent shock to investment specific technical change"
Eps_mu "Transitory shock to investment specific technical change"
Eps_xitheta "Collateral shock"

parameters
%--- 5 estimated parameters
gamma_h gamma_e omega g_trans lambda_qbar_trans
%--- 8 calibrated parameters
alpha r ky ik qlLeOY theta_bar qlLhOY n lbar
%--- 8 parameters related to shock processes
rho_z rho_v rho_q rho_mu rho_a rho_phi rho_psi rho_xitheta

sig_Eps_a, sig_Eps_psi, sig_Eps_z, sig_Eps_v, sig_Eps_q, sig_Eps_mu, sig_Eps_xitheta

hetero_tp_1_2, hetero_tp_2_1

parameters(hetero,2) sig_Eps_phi


%--------------------------------------------------------------------
%-- Model specification
%--------------------------------------------------------------------
model
	# g = g_trans/100 + 1;
	# lambda_qbar = 0.01*lambda_qbar_trans + 1;
	# lambda_kbar = g*lambda_qbar; %See second sentence of first full paragraph of page 9
	# delta = 1-(1-ik)*lambda_kbar ;  %Rearrange equation 60 
	
	%From Tao's calibration5.m file
	# betaparam = (qlLeOY*(1.0-g*theta_bar/r)+ky*(1.0-g*theta_bar/(lambda_kbar*r)))/(alpha+qlLeOY*(1.0-theta_bar)+ky*(1.0-delta-theta_bar)/lambda_kbar); 
	# lambda_abar = g/(betaparam*r)-1.0; %See equation 58 (first part)
	# muboe = (betaparam*lambda_abar)/g; %See equation 58 (second part)
	# phi = qlLeOY*(1.0-betaparam*(1.0+lambda_abar*theta_bar))/(betaparam*alpha); % from calibration5.m
	% # qlLeOY = (betaparam*alpha*phi)/(1 - betaparam - betaparam*lambda_abar*theta_bar) %Equation 59;
	%Do not need equation 61
	
	# iy = ik*ky; %Equation 62
	# by = theta_bar*g*qlLeOY + (theta_bar/lambda_qbar)*ky; %Equation 63
	
	
	# ceOY = alpha-iy-by*(1-betaparam*(1+lambda_abar))/g; %See equation 64
	# chOY = 1-ceOY-iy; %See equation 65
	
	# lhOle = qlLhOY/qlLeOY; %See equation 67
	# le = lbar/(1.0+lhOle);
	# lh = lbar-le;
	
	# psi_bar = ((1-alpha)*g*(1-gamma_h/r)*(1/chOY))/(n*(g-gamma_h)); %Equation 68 
	
	# omegah = (g - betaparam*(1+lambda_abar)*gamma_h)*(g - gamma_h); %See top of page 12
	# omegae = (g - betaparam*gamma_e)*(g - gamma_e); %See top of page 12
	
	%(1)
	omegah*Muh = -(g^2 + (gamma_h^2)*betaparam*(1+lambda_abar))*Ch + g*gamma_h*(Ch(-1) - G_Lambda) -
	betaparam*lambda_abar*gamma_h*(g - gamma_h)*Lambda_a(+1) + 
	betaparam*(1+lambda_abar)*g*gamma_h*(Ch(+1) + G_Lambda(+1)); %Eqn (69)
	
	
	%(2)
	W + Muh = Psi ; % Eqn. 70
	
	%(3)
	Ql + Muh = betaparam*(1+lambda_abar)*(Muh(+1) + Ql(+1)) + (1 - betaparam*(1+lambda_abar))*(Phi -
	Lh) + betaparam*lambda_abar*Lambda_a(+1); % Eqn. 71
	
	
	%(4) 
	Muh - R = Muh(+1) + (lambda_abar/(1+lambda_abar))*Lambda_a(+1) - G_Lambda(+1); % Eqn (72)
	
	%(5)
	omegae*Mue = -(g^2 + betaparam*gamma_e^2)*Ce + g*gamma_e*(Ce(-1) - G_Lambda) +
	betaparam*g*gamma_e*(Ce(+1) + G_Lambda(+1)); % Eqn (73)
	
	%(6)
	W = Y - N; % Eqn (74)
	
	%(7)
	Qk = (1+betaparam)*omega*(lambda_kbar^2)*I - omega*(lambda_kbar^2)*I(-1) +
	omega*(lambda_kbar^2)*(G_Lambda + G_Q) - betaparam*omega*(lambda_kbar^2)*(I(+1)
	 + G_Lambda(+1) + G_Q(+1)); % Eqn (75)
	
	%(8)
	Qk + Mue = (muboe*theta_bar/lambda_qbar)*(Mub + Thetat) + 
	(betaparam*(1-delta)/(lambda_kbar))*(Qk(+1) - G_Q(+1) - G_Lambda(+1)) + 
	(1 - (muboe*theta_bar/lambda_qbar))*Mue(+1) + 
	(muboe*theta_bar/lambda_qbar)*(Qk(+1) - G_Q(+1)) + 
	betaparam*alpha*(1 - phi)*(1/ky)*(Y(+1) - K); % Eqn (76)
	
	%(9)
	Ql + Mue = muboe*g*theta_bar*(Thetat + Mub) + (1 - muboe*g*theta_bar)*Mue(+1) + 
	muboe*g*theta_bar*(Ql(+1) + G_Lambda(+1)) + betaparam*Ql(+1) + 
	(1 - betaparam - betaparam*lambda_abar*theta_bar)*(Y(+1) - Le); % Eqn (77)
	
	%(10)
	Mue - R = (1/(1+lambda_abar))*(Mue(+1) - G_Lambda(+1) + lambda_abar*Mub); % Eqn (78)
	
	%(11)
	Y = alpha*phi*Le(-1) + alpha*(1- phi)*K(-1) + (1-alpha)*N - 
	(((1-phi)*alpha)/(1-((1-phi)*alpha)))*(G_Z + G_Q); % Eqn (79)
	
	%(12)
	K = ((1-delta)/(lambda_kbar))*(K(-1) - G_Lambda - G_Q) 
	+ (1 - ((1-delta)/(lambda_kbar)))*I; % Eqn (80)
	
	%(13)
	Y = chOY*Ch + ceOY*Ce + iy*I; % Eqn (81)
	
	%(14)
	0 = (lh/lbar)*Lh + (le/lbar)*Le;  % Eqn (82)
	
	%(15)
	alpha*Y = ceOY*Ce + iy*I + qlLeOY*(Le - Le(-1)) + ((1/g)*by)*(B(-1) - G_Lambda) 
	- ((1/r)*by)*(B - R);% Eqn (83)
	
	%(16)
	B = Thetat + (g*theta_bar*(qlLeOY/by))*(Ql(+1) + Le + G_Lambda(+1)) + 
	(1 - (g*theta_bar*(qlLeOY/by)))*(Qk(+1) + K - G_Q(+1)); % Eqn (84)
	
	%(17)
	G_Z = Lambda_z + V - V(-1); % Eqn (85)
	%(18)
	G_Q = Lambda_q + Mu - Mu(-1); % Eqn (86)
	%(19)
	G_Lambda = (1/(1-(1-phi)*alpha))*G_Z + (((1-phi)*alpha)/(1-(1-phi)*alpha))*G_Q; % Eqn (87)
	
	
	%Shock processes
	
	Lambda_z = rho_z*Lambda_z(-1) + sig_Eps_z*Eps_z; % Eqn (88)
	
	V = rho_v*V(-1) + sig_Eps_v*Eps_v; % Eqn (89)
	
	Lambda_q = rho_q*Lambda_q(-1) + sig_Eps_q*Eps_q; % Eqn (90)
	
	Mu = rho_mu*Mu(-1) + sig_Eps_mu*Eps_mu; % Eqn (91)
	
	Lambda_a = rho_a*Lambda_a(-1) + sig_Eps_a*Eps_a; % Eqn (93)
	
	Phi = rho_phi*Phi(-1) + sig_Eps_phi*Eps_phi; % Eqn (94)
	
	Psi = rho_psi*Psi(-1) + sig_Eps_psi*Eps_psi; % Eqn (95)
	
	Thetat = rho_xitheta*Thetat(-1) + sig_Eps_xitheta*Eps_xitheta; %Eqn (96)
	
	C = (chOY/(chOY+ceOY))*Ch + (ceOY/(chOY+ceOY))*Ce;
	
	%--- Observer equations
	DLogQl = Ql - Ql(-1) + log(g) + G_Lambda;
	DLogQ =  G_Q + log(lambda_qbar);
	DLogC = (chOY/(chOY+ceOY))*(Ch - Ch(-1)) + (ceOY/(chOY+ceOY))*(Ce - Ce(-1)) + log(g) + G_Lambda; 
	DLogI = I - I(-1) + log(g) + G_Lambda; 
	DLogB = B - B(-1) + log(g) + G_Lambda;
	LogL = N + log(n);
	
	%--- Auxillary variables
	DLogGDP = iy*DLogI + (1-iy)*DLogC;

steady_state_model	
	DLogQl = log(g);
	DLogQ = log(lambda_qbar);
	DLogC = log(g); 
	DLogI = log(g); 
	DLogB = log(g);
	LogL = log(n);
	DLogGDP = iy*DLogI + (1-iy)*DLogC;
	