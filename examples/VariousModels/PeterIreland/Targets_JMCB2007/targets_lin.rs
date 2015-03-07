% Peter N. Ireland (2007): "Changes in the Federal Reserve's Inflation Target:
% Causes and Consequences" Journal of Money, Credit and Banking, Vol. 39, No. 8 (December 2007)
%
% linear model
% data start 1959Q1

endogenous Y "Output" PAI "Gross inflation" PAISTAR "Inflation target" A LAMBDA Z
R "Gross nom. int. rate" GY GPAI GR RRPAI Q X "Output gap" V E
GY_HAT GPAI_HAT RRPAI_HAT

observables GY_HAT GPAI_HAT RRPAI_HAT

exogenous EPS_A "Preference" EPS_E "Cost-push" EPS_Z "Technology" EPS_PAI "Inflation target" EPS_V "Monetary policy"

parameters	psi	alpha rho_a	sig_a gam beta rho_e sig_e zss sig_z rho_pai rho_x rho_gy
	delta_a	delta_e	delta_z	sig_pai	rho_v sig_v

model
	
	 A = rho_a*A{-1}+sig_a*EPS_A;

	 (zss-gam)*(zss-beta*gam)*LAMBDA = gam*zss*Y{-1}-(zss^2+beta*gam^2)*Y+
	 	beta*gam*zss*Y{+1}+(zss-gam)*(zss-beta*gam*rho_a)*A-gam*zss*Z;

	 LAMBDA = LAMBDA{+1}+R-PAI{+1};

	 E = rho_e*E{-1}+sig_e*EPS_E;

	 Z=sig_z*EPS_Z;

	 (1+beta*alpha)*PAI = alpha*PAI{-1}+beta*PAI{+1}-psi*LAMBDA+psi*A-E-alpha*PAISTAR;

	 GY = Y-Y{-1}+Z;

	 GPAI = PAI-PAI{-1}+PAISTAR;

	 GR = R-R{-1}+PAISTAR;

	 RRPAI = R-PAI;

	 0 = gam*zss*Q{-1} -(zss^2+beta*gam^2)*Q+beta*gam*zss*Q{+1}+beta*gam*(zss-gam)*(1-rho_a)*A-gam*zss*Z;

	 X = Y-Q;

	 R = R{-1}+rho_pai*PAI+rho_x*X+rho_gy*GY-PAISTAR+V;

	 PAISTAR=delta_a*EPS_A-delta_e*EPS_E-delta_z*EPS_Z+sig_pai*EPS_PAI;

	 V=rho_v*V{-1}+sig_v*EPS_V;

	 % Measurement equations
	 %----------------------

	 GY_HAT=GY;

	 GPAI_HAT=GPAI;

	 RRPAI_HAT=RRPAI;
	