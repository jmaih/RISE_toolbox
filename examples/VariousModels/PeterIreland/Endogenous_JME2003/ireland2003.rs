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

	g*K = (1-delta)*K+X*I;

	log(X) = rho_x*log(X{-1})+sig_x*EPS_X;

	Y = C+I+phi_k/2*(K/K{-1}-1)^2*K{-1}+phi_p/2*(PAI/steady_state(PAI)-1)^2*Y;

	log(A) = rho_a*log(A{-1})+sig_a*EPS_A;

	log(E) = (1-rho_e)*log(e_ss)+rho_e*log(E{-1})+sig_e*EPS_E;

	A = LAMBDA*C^(1/gam)*(C^((gam-1)/gam)+E^(1/gam)*M^((gam-1)/gam));

	eta = LAMBDA*W*(1-H);

	E*C = M*(1-1/R)^gam;

	g*LAMBDA =beta*R*LAMBDA{+1}/PAI{+1};

	g*LAMBDA/X+phi_k*LAMBDA*(K/K{-1}-1) = beta*LAMBDA{+1}*(Q{+1}+(1-delta)/X{+1})
		-beta*phi_k/2*LAMBDA{+1}*(K{+1}/K-1)^2
		+beta*phi_k*LAMBDA{+1}*(K{+1}/K-1)*K{+1}/K;
		
	log(Z) = (1-rho_z)*log(z_ss)+rho_z*log(Z{-1})+sig_z*EPS_Z;

	D = Y-W*H-Q*K-phi_p/2*(PAI/steady_state(PAI)-1)^2*Y;

	Y = K{-1}^alpha*(Z*H)^(1-alpha);

	LAMBDA*W*H = (1-alpha)*XI*Y;

	LAMBDA*Q*K = alpha*XI*Y;

	phi_p*LAMBDA*(PAI/steady_state(PAI)-1)*PAI/steady_state(PAI) = (1-theta)*LAMBDA+theta*XI+
		beta*phi_p*LAMBDA{+1}*(PAI{+1}/steady_state(PAI)-1)*PAI{+1}/steady_state(PAI)*Y{+1}/Y;

	omega_r*log(R/steady_state(R)) = omega_mu*log(MU/mu_ss)+omega_pai*log(PAI/steady_state(PAI))+
		omega_y*log(Y/steady_state(Y))+log(V);

	log(V) = rho_v*log(V{-1})+sig_v*EPS_V;

	N = Y/H;

	MU =g*M/M{-1}*PAI;

	% Measurement equations
	%-----------------------

	LC=log(C)+log(g)*TREND;

	LI=log(I)+log(g)*TREND;

	LM=log(M)+log(g)*TREND;

	LPI=log(PAI);

	LR=log(R);

steady_state_model	 % Y C I H N M MU
	X=1;
	A=1;
	E=e_ss;
	Z=z_ss;
	V=1;
	% mu determined by policy
	MU=mu_ss;
	PAI=MU/g;
	R=MU/beta;
	Q=g/beta-1+delta;
	% lambda comes here
	LAMBDA=(
		eta+(1-alpha)*1/(1+E*(R/(R-1))^(gam-1))*1/(theta/(theta-1)-(g-1+delta)*alpha/Q)
		)
		/(
		(1-alpha)*Z*((theta-1)/theta)^(1/(1-alpha))*(alpha/Q)^(alpha/(1-alpha))
		);
	XI = (theta-1)/theta*LAMBDA;
	C = 1/(1+E*(R/(R-1))^(gam-1))*(1/LAMBDA);
	M = E*(R/(R-1))^gam*C;
	Y = 1/(1-(g-1+delta)*(theta-1)/theta*alpha/Q)*C;
	K = (theta-1)/theta*alpha*Y/Q;
	I = (g-1+delta)*K;
	H=1/Z*(Y/K^alpha)^(1/(1-alpha));
	W=(1-alpha)*(theta-1)/theta*Y/H;
	D = Y-W*H-Q*K;
	N=Y/H;
	LC=log(C);
	LI=log(I);
	LM=log(M);
	LPI=log(PAI);
	LR=log(R);

parameterization
	% fixed parameters
	%------------------
	g		  ,1 ;
	delta	  ,0.025 ;
	eta		  ,1.5 ;
	theta	  ,6 ;
	sig_v	  ,0.01 ;
	% estimated parameters
	%-----------------------
	beta	   ,0.9980 , 0.9, 1;
	gam		   ,0.0736 , 0.0005, 2 ;
	alpha	   ,0.2022 , 0.0005, 1 ;
	phi_p_trans,54.0745/100, 0, 2;
	phi_k_trans,12.4368/100, 0, 2;
	mu_ss_trans,1.0110-1, 0.0005, 2 ;
	omega_r	   ,3.0296 , 0.5, 5;
	omega_mu   ,0.9840 , -3,  3;
	omega_y	   ,-0.0239 , -3,  3 ;
	omega_pai  ,2.0070 , -3,  3 ;
	e_ss	   ,2.7599 , 0.0005, 5  ;
	z_ss_trans ,7034.6/10000, 0, 2;
	rho_a	   ,0.9903 , 0, 1 ;
	rho_e	   ,0.9497 , 0, 1 ;
	rho_x	   ,0.6975 , 0, 1 ;
	rho_z	   ,0.9787 , 0, 1 ;
	rho_v	   ,0.4400 , 0, 1 ;
	sig_a	   ,0.0064 , 0.0005, 2 ;
	sig_e	   ,0.0115 , 0.0005, 2  ;
	sig_x	   ,0.0224 , 0.0005, 2  ;
	sig_z	   ,0.0153 , 0.0005, 2  ;
	