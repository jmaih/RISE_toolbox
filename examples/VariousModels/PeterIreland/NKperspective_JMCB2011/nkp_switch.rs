% Peter N. Ireland (2011): "A New Keynesian Perspective on the Great Recession"
% Journal of Money, Credit and Banking, Vol 43, No. 1
%
% Nonlinear model
% Data starting in 1983:1

endogenous Y "Output" C "Consumption" PAI "Inflation" A LAMBDA Z THETA G "Output growth" Q
X "Output gap" R "interest rate" E GHAT PAIHAT RHAT

observables GHAT PAIHAT RHAT

exogenous EPS_A "Preference" EPS_Z "Technology" EPS_R "monetary policy" EPS_E "Cost push"

parameters psi alpha rho_a beta gam rho_e thetass zss paiss rho_r ess  pol_tp_1_2  pol_tp_2_1  vol_tp_1_2  vol_tp_2_1

parameters(pol,2) rho_pai rho_x	rho_g

parameters(vol, 2) sig_a sig_e sig_z sig_r

model

	# rho_theta = rho_e;
	# phi = (thetass-1)/psi;
	# sig_theta = sig_e*phi;

	Y = C + phi/2*(PAI/(PAI{-1}^alpha*paiss^(1-alpha))-1)^2*Y;

	log(A) = rho_a*log(A{-1})+sig_a*EPS_A;

	LAMBDA = A*Z/(Z*C-gam*C{-1})-beta*gam*(A{+1}/(Z{+1}*C{+1}-gam*C));

	LAMBDA = beta*R*LAMBDA{+1}/(Z{+1}*PAI{+1});

	log(E/ess) = rho_e*log(E{-1}/ess)+sig_e*EPS_E;

%	log(THETA) = (1-rho_theta)*log(thetass)+rho_theta*log(THETA{-1})+sig_theta*EPS_THETA;

	log(Z) = log(zss) + sig_z*EPS_Z;

	THETA-1 = THETA*A/LAMBDA-phi*(PAI/(PAI{-1}^alpha*paiss^(1-alpha))-1)*(PAI/(PAI{-1}^alpha*paiss^(1-alpha)))
			  +beta*phi*LAMBDA{+1}*Y{+1}/(LAMBDA*Y)*(PAI{+1}/(PAI^alpha*paiss^(1-alpha))-1)*
			  PAI{+1}/(PAI^alpha*paiss^(1-alpha));

	G = Y/Y{-1}*Z;

	1 = Z/(Z*Q-gam*Q{-1})-beta*gam*A{+1}/A*1/(Z{+1}*Q{+1}-gam*Q);

	X = Y/Q;

	log(R/steady_state(R)) =rho_r*log(R{-1}/steady_state(R))+rho_pai*log(PAI/paiss)+
		rho_x*log(X/steady_state(X))+rho_g*log(G/steady_state(G))+sig_r*EPS_R;

	THETA=thetass*(E/ess)^(-phi);

	% Measurement equations
	%----------------------
	GHAT = log(G/steady_state(G));
	
	PAIHAT = log(PAI/steady_state(PAI));
	
	RHAT = log(R/steady_state(R));

steady_state_model

	A=1;

	THETA=thetass;

	Z = zss;

	PAI = paiss;

	LAMBDA = thetass/(thetass-1);

	Y =(thetass-1)/thetass*(zss-beta*gam)/(zss-gam);

	C = Y;

	G = Z;

	Q = (zss-beta*gam)/(zss-gam);

	X= (thetass-1)/thetass;

	R=zss/beta*PAI;

	E = ess;
