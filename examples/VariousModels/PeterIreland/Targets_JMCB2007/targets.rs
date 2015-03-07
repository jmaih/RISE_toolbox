% Peter N. Ireland (2007): "Changes in the Federal Reserve's Inflation Target:
% Causes and Consequences" Journal of Money, Credit and Banking, Vol. 39, No. 8 (December 2007)
%
% Nonlinear model
% data start 1959Q1

endogenous Y "Output" C "Consumption" PAI "Gross inflation" PAISTAR "Inflation target" A LAMBDA Z
R "Gross nom. int. rate" THETA GY GPAI GR RRPAI Q X "Output gap" V E EPS_THETA
GY_HAT GPAI_HAT RRPAI_HAT

observables GY_HAT GPAI_HAT RRPAI_HAT

exogenous EPS_A "Preference" EPS_E "Cost-push" EPS_Z "Technology" EPS_PAI "Inflation target" EPS_V "Monetary policy"

parameters	psi	alpha rho_a	sig_a gam beta rho_e thetass sig_e zss sig_z rho_pai rho_x rho_gy
	delta_a	delta_e	delta_z	sig_pai	rho_v sig_v ess

model
	% Auxiliary parameters (see definitions on page 9 in the notes)
	%--------------------------------------------------------------
	# phi =(thetass-1)/psi;

	# sig_theta=phi*sig_e;
	
	# rho_theta=rho_e;

	# delta_theta = delta_e;
	
	 Y = C +phi/2*(PAI*(PAISTAR/PAI{-1})^alpha-1)^2*Y;

	 log(A) = rho_a*log(A{-1})+sig_a*EPS_A;

	 LAMBDA = A*Z/(Z*C-gam*C{-1})-beta*gam*A{+1}/(Z{+1}*C{+1}-gam*C);

	 LAMBDA = beta*R*1/Z{+1}*1/PAISTAR{+1}*LAMBDA{+1}/PAI{+1};

	 log(THETA) = (1-rho_theta)*log(thetass)+rho_theta*log(THETA{-1})+sig_theta*EPS_THETA;

	 log(Z)=log(zss)+sig_z*EPS_Z;

	 THETA-1 = THETA*A/LAMBDA-phi*(PAI*(PAISTAR/PAI{-1})^alpha-1)*PAI*(PAISTAR/PAI{-1})^alpha+
	 		beta*phi*LAMBDA{+1}/LAMBDA*(PAI{+1}*(PAISTAR{+1}/PAI)^alpha-1)*PAI{+1}*(PAISTAR{+1}/PAI)^alpha*Y{+1}/Y;

	GY = Y/Y{-1}*Z;

	GPAI=PAI/PAI{-1}*PAISTAR;

	GR=R/R{-1}*PAISTAR;

	RRPAI = R/PAI;

	1 = Z/(Z*Q-gam*Q{-1})-beta*gam*A{+1}/A*1/(Z{+1}*Q{+1}-gam*Q);

	X = Y/Q;

	log(R) = log(R{-1})+rho_pai*log(PAI)+rho_x*log(X/steady_state(X))+
		rho_gy*log(GY/steady_state(GY))-log(PAISTAR)+log(V);

	log(PAISTAR)=delta_a*EPS_A-delta_theta*EPS_THETA-delta_z*EPS_Z+sig_pai*EPS_PAI;

	log(V)=rho_v*log(V{-1})+sig_v*EPS_V;

	% Auxiliary equations (page 9 in the notes)
	%-------------------------------------------

	log(E/ess) = 1/phi*log(THETA/thetass);

	EPS_THETA=EPS_E;

	% Measurement equations
	%----------------------

	GY_HAT=log(GY/steady_state(GY));

	GPAI_HAT=log(GPAI/steady_state(GPAI));

	RRPAI_HAT=log(RRPAI/steady_state(RRPAI));

%steady_state_model
%
%	A=1;
%	THETA=thetass;
%	Z=zss;
%	PAISTAR=1;
%	V=1;
%	PAI=1;
%	Y=(thetass-1)/thetass*(zss-beta*gam)/(zss-gam);
%	C=Y;
%	GY=Z;
%	GPAI=1;
%	GR=1;
%	LAMBDA=thetass/(thetass-1);
%	Q=(zss-beta*gam)/(zss-gam);
%	X=(thetass-1)/thetass;
%	R=zss/beta;
%	RRPAI=zss/beta;
%	E=ess;
