%
%This file replicates the estimation of the cash in advance model described
%Frank Schorfheide (2000): "Loss function-based evaluation of DSGE models",
%Journal of Applied Econometrics, 15(6), 645-670.
%
%The data are in file "fsdat_simul.m", and have been artificially generated.
%They are therefore different from the original dataset used by Schorfheide.
%
%The equations are taken from J. Nason and T. Cogley (1994): "Testing the
%implications of long-run neutrality for monetary business cycle models",
%Journal of Applied Econometrics, 9, S37-S70.
%Note that there is an initial minus sign missing in equation (A1), p. S63.
%
%This implementation was written by Michel Juillard. Please note that the
%following copyright notice only applies to this Dynare implementation of the
%model.
%

endogenous m P c e W R k d n l gy_obs gp_obs y dA

exogenous e_a e_m

parameters alp "$\alpha$" bet "$\beta$" gam "$\gamma$" mst rho "$\rho$" psi "$\psi$"
	del "$\delta$" sig_a "$\sigma_a$" sig_m  "$\sigma_m$"

observables gp_obs gy_obs

model
	dA = exp(gam+sig_a*e_a);
	
	log(m) = (1-rho)*log(mst) + rho*log(m(-1))+sig_m*e_m;
	
	-P/(c(+1)*P(+1)*m)+bet*P(+1)*(alp*exp(-alp*(gam+log(e(+1))))*k^(alp-1)*n(+1)^(1-alp)+(1-del)*exp(-(gam+log(e(+1)))))/(c(+2)*P(+2)*m(+1))=0;
	
	W = l/n;
	
	-(psi/(1-psi))*(c*P/(1-n))+l/n = 0;
	
	R = P*(1-alp)*exp(-alp*(gam+sig_a*e_a))*k(-1)^alp*n^(-alp)/W;
	
	1/(c*P)-bet*P*(1-alp)*exp(-alp*(gam+sig_a*e_a))*k(-1)^alp*n^(1-alp)/(m*l*c(+1)*P(+1)) = 0;
	
	c+k = exp(-alp*(gam+sig_a*e_a))*k(-1)^alp*n^(1-alp)+(1-del)*exp(-(gam+sig_a*e_a))*k(-1);
	
	P*c = m;
	
	m-1+d = l;
	
	e = exp(sig_a*e_a);
	
	y = k(-1)^alp*n^(1-alp)*exp(-alp*(gam+sig_a*e_a));
	

% We have to take into account the loglinear option of dynare
% (which will take the log of the steady state during estimation)	
% In rise, we have to explicitly put the variables in log terms
	exp(gy_obs) = dA*y/y(-1);  %	gy_obs = dA*y/y(-1);
	
	exp(gp_obs) = (P/P(-1))*m(-1)/dA;  %	gp_obs = (P/P(-1))*m(-1)/dA;

	
%steady_state_model;	% if this does not work, rise will use it as initial guess
%	k = 6;
%	m = mst;
%	P = 2.25;
%	c = 0.45;
%	e = 1;
%	W = 4;
%	R = 1.02;
%	d = 0.85;
%	n = 0.19;
%	l = 0.86;
%	y = 0.6;
%	gy_obs = gam;  %	gy_obs = exp(gam);
%	gp_obs = -gam; %	gp_obs = exp(-gam);
%	dA = exp(gam);
%end;

%xx_ssmdef_1 ... xx_ssmdef_9 are known words to rise
steady_state_model
  dA = exp(gam);
  xx_ssmdef_1 = 1/dA;
  m = mst;
  xx_ssmdef_2 = ( (1-xx_ssmdef_1*bet*(1-del)) / (alp*xx_ssmdef_1^alp*bet) )^(1/(alp-1));
  xx_ssmdef_3 = ( ((xx_ssmdef_2*xx_ssmdef_1)^alp - (1-xx_ssmdef_1*(1-del))*xx_ssmdef_2)/mst )^(-1);
  xx_ssmdef_4 = psi*mst^2/( (1-alp)*(1-psi)*bet*xx_ssmdef_1^alp*xx_ssmdef_2^alp );
  n  = xx_ssmdef_3/(xx_ssmdef_4+xx_ssmdef_3);
  P  = xx_ssmdef_3 + xx_ssmdef_4;
  k  = xx_ssmdef_2*n;

  l  = psi*mst*n/( (1-psi)*(1-n) );
  c  = mst/P;
  d  = l - mst + 1;
  y  = k^alp*n^(1-alp)*xx_ssmdef_1^alp;
  R  = mst/bet;
  W  = l/n;
%%  ist  = y-c;
%  q  = 1 - d;

  e = 1;
  
  gp_obs = log(m/dA);% accommodating dynare's loglinear option
  gy_obs = log(dA); % accommodating dynare's loglinear option

parameterization; % point at which dynare starts estimation 								% dynare calibration
	alp  ,   0.356000,   0.3235,    0.3891,  beta_pdf(.90)     ;% 0.356, 0.02;	      	0.330
	bet  ,   0.993000,   0.9896,    0.9958,  beta_pdf(.90)     ;% 0.993, 0.002;			0.990
	gam  ,   0.008500,   0.0036,    0.0134,  normal_pdf(.90)   ;% 0.0085, 0.003;		  	0.003
	mst  ,   1.000200,   0.9887,    1.0117,  normal_pdf(.90)   ;% 1.0002, 0.007;		  	1.011
	rho  ,   0.129000,   0.0001,    0.6851,  beta_pdf(.90)     ;% 0.129, 0.223;			0.700
	psi  ,   0.650000,   0.5658,    0.7304,  beta_pdf(.90)     ;% 0.65, 0.05;			0.787
	del  ,   0.010000,   0.0034,    0.0194,  beta_pdf(.90)     ;% 0.01, 0.005;			0.020
	sig_a,   0.035449,   0.0075,    0.0998,  inv_gamma_pdf(.90);% 0.035449, inf;		  	0.014
	sig_m,   0.008862,   0.0019,    0.0249,  inv_gamma_pdf(.90);% 0.008862, inf;		  	0.005

% Look into folder Moments2Bounds to see how to go from moments to bounds of the distribution
