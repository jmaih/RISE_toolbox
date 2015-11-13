endogenous Y C I H N M W Q D P R K LAMB XI A E Z X V YGT PAI MU

log_vars Y C I H N M W Q D P R K LAMB XI A E Z X V YGT PAI MU

exogenous EA EE TREND EZ EV EX

parameters delta rho_x sig_x phi_k_trans g phi_p_trans rho_a sig_a rho_e e_ss sig_e gam eta beta rho_z sig_z

	omega_r omega_mu omega_pai omega_y rho_v sig_v z_ss_trans theta alpha mu_ss_trans

model

	# z_ss=z_ss_trans*10000;
	# phi_p = 100*abs(phi_p_trans);
	# phi_k = 100*abs(phi_k_trans);
	# mu_ss = 1 + mu_ss_trans;

	YGT = Y/g^TREND;

	PAI = P/P{-1};

	MU = M/M{-1};

	K = (1-delta)*K{-1}+X*I;

	log(X) = rho_x*log(X{-1})+sig_x*EX;

	Y = C + I + phi_k/2*(K/(g*K{-1})-1)^2*K{-1} + phi_p/2*(P/(steady_state(PAI)*P{-1})-1)^2*Y;

	log(A) = rho_a*log(A{-1})+sig_a*EA;

	log(E) = (1-rho_e)*log(e_ss)+rho_e*log(E{-1})+sig_e*EE;

	A = LAMB*C^(1/gam)*(C^((gam-1)/gam)+E^(1/gam)*(M/P)^((gam-1)/gam));

	eta = LAMB*W/P*(1-H);

	E*C = M/P*(1-1/R)^gam;

	LAMB = beta*R*LAMB{+1}*P/P{+1};

	LAMB*(1/X+phi_k/g*(K/(g*K{-1})-1)) = beta*LAMB{+1}*(Q{+1}/P{+1}+(1-delta)/X{+1}) -
		beta*phi_k/2*LAMB{+1}*(K{+1}/(g*K)-1)^2 +
		beta*phi_k*LAMB{+1}*(K{+1}/(g*K)-1)*K{+1}/(g*K);

	log(Z) = (1-rho_z)*log(z_ss)+rho_z*log(Z{-1})+sig_z*EZ;

	D/P = Y -(W*H+Q*K{-1})/P-phi_p/2*(P/(steady_state(PAI)*P{-1})-1)^2*Y;

	Y = K{-1}^alpha*(g^TREND*Z*H)^(1-alpha);

	LAMB*W/P*H = (1-alpha)^XI*Y;

	LAMB*Q/P*K{-1} = alpha*XI*Y;

	phi_p*LAMB*(P/(steady_state(PAI)*P{-1})-1)*P/(steady_state(PAI)*P{-1}) = (1-theta)*LAMB+ theta*XI +
		beta*phi_p*LAMB{+1}*(P{+1}/(steady_state(PAI)*P)-1)*P{+1}/(steady_state(PAI)*P)*Y{+1}/Y;

	omega_r*log(R/steady_state(R)) = omega_mu*log(MU/mu_ss) + omega_pai*log(PAI/steady_state(PAI)) +
		omega_y*log(YGT/steady_state(YGT))+log(V);

	log(V) = rho_v*log(V{-1}) + sig_v*EV;

	N = Y/H;

	