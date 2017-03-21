% Peter N. Ireland (2003): "Endogenous money of sticky prices?"
% Journal of Monetary Economics 50 (2003) 1623--1648
%
% Log-linearized version of the model with a flag "original" for choosing between
% the original equations or the slightly rewritten ones. (see page 13 in the notes)
%
% The alternative form allows Peter Ireland to get rid of variables that
% are both predetermined and forward looking. But we don't need to do this
% in RISE

endogenous Y C I H N M MU W Q D PAI R K LAMBDA XI A E Z X V
LC LI LM LPI LR

observables LC LI LM LPI LR

exogenous EPS_V EPS_Z EPS_E EPS_A EPS_X	TREND

parameters g delta rho_x sig_x phi_p_trans phi_k_trans
 rho_a sig_a rho_e e_ss sig_e z_ss_trans
		   gam eta beta	rho_z  sig_z alpha theta mu_ss_trans rho_v sig_v
		   omega_r omega_mu	omega_pai omega_y
		   
model

	# z_ss=z_ss_trans*10000;
	# phi_p = 100*abs(phi_p_trans);
	# phi_k = 100*abs(phi_k_trans);
	# mu_ss = 1 + mu_ss_trans;
	# xss=1;
    # ass=1;
    # ess=e_ss;
    # zss=z_ss;
    # vss=1;
    # muss=mu_ss;
    # paiss=muss/g;
    # rss=muss/beta;
    # qss=g/beta-1+delta;
    # l1ss=1/(theta/(theta-1)-(g-1+delta)*alpha/qss);
    # l2ss=1/(1+ess*(rss/(rss-1))^(gam-1));
    # l3ss=(1-alpha)*zss*((theta-1)/theta)^(1/(1-alpha))*(alpha/qss)^(alpha/(1-alpha));
    # lambdass=(eta+(1-alpha)*l1ss*l2ss)/l3ss;
    # xiss = (theta-1)/theta*lambdass;
    # css = 1/(1+ess*(rss/(rss-1))^(gam-1))*(1/lambdass);
    # mss = ess*(rss/(rss-1))^gam*css;
    # yss = 1/(1-(g-1+delta)*(theta-1)/theta*alpha/qss)*css;
    # kss = (theta-1)/theta*alpha*yss/qss;
    # iss = (g-1+delta)*kss;
    # hss=1/zss*(yss/kss^alpha)^(1/(1-alpha));
    # wss=(1-alpha)*(theta-1)/theta*yss/hss;
    # dss = yss-wss*hss-qss*kss;
    # nss=yss/hss;
    # lcss=log(css);
    # liss=log(iss);
    # lmss=log(mss);
    # lpiss=log(paiss);
    # lrss=log(rss);

	%-----------------------------------------------------
	g*kss*K = (1-delta)*kss*K{-1}+iss*X+iss*I;

	X = rho_x*X{-1}+sig_x*EPS_X;

	yss*Y = css*C+iss*I;

	A = rho_a*A{-1}+sig_a*EPS_A;

	E = rho_e*E{-1}+sig_e*EPS_E;

	@#if original
		gam*rss*A = gam*rss*LAMBDA+rss*(1+(gam-1)*lambdass*css)*C+(rss-1)*lambdass*mss*E+
					(gam-1)*(rss-1)*lambdass*mss*M;
	@#else
		gam*rss*A = gam*rss*LAMBDA+rss*(1+(gam-1)*lambdass*css)*C+(rss-1)*lambdass*mss*E
					+(gam-1)*(rss-1)*lambdass*mss*MU
					+(gam-1)*(rss-1)*lambdass*mss*M{-1}
					-(gam-1)*(rss-1)*lambdass*mss*PAI;
	@#end

	0 = eta*LAMBDA+eta*W-lambdass*wss*hss*H;

	@#if original
		(rss-1)*E+(rss-1)*C = (rss-1)*M+gam*R;
	@#else
		(rss-1)*E+(rss-1)*C = (rss-1)*MU+(rss-1)*M{-1}-(rss-1)*PAI+gam*R;
	@#end

	LAMBDA =R+LAMBDA{+1}-PAI{+1};

	@#if original
		g*LAMBDA-g*X-phi_k*K{-1} = g*LAMBDA{+1}+beta*qss*Q{+1}-beta*(1-delta)*X{+1}
			+beta*phi_k*K{+1}
			-(1+beta)*phi_k*K;
	@#else
		g*LAMBDA-(g+beta*(phi_k*(1-(1-delta)/g)-(1-delta))*rho_x)*X-phi_k*K{-1} = 
		g*LAMBDA{+1}+beta*qss*Q{+1}+(beta*(1-delta)/g-(1+beta))*phi_k*K
		+beta*phi_k*(1-(1-delta)/g)*I{+1};
	@#end
		
	Z = rho_z*Z{-1}+sig_z*EPS_Z;

	dss*D = yss*Y-wss*hss*W-wss*hss*H-qss*kss*Q-qss*kss*K{-1};

	Y = alpha*K{-1}+(1-alpha)*Z+(1-alpha)*H;

	LAMBDA+W+H = XI+Y;

	LAMBDA+Q+K{-1} = XI+Y;

	phi_p*PAI = (1-theta)*LAMBDA+(theta-1)*XI+beta*phi_p*PAI{+1};

	omega_r*R = omega_mu*MU+omega_pai*PAI+omega_y*Y+V;

	V = rho_v*V{-1}+sig_v*EPS_V;

	N = Y-H;

	MU =M-M{-1}+PAI;

	% Measurement equations
	%-----------------------

	LC=C+log(css)+log(g)*TREND;

	LI=I+log(iss)+log(g)*TREND;

	LM=M+log(mss)+log(g)*TREND;

	LPI=PAI+log(paiss);

	LR=R+log(rss);