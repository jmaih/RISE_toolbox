%-------------------------------------------------------------
% Nonlinear stationary RBC model
% Reference: Foerster, Rubio-Ramirez, Waggoner and Zha (2014)
% Perturbation Methods for Markov-Switching Dynamic Stochastic
% General Equilibrium Models.
% Quantitative Economics 7(2016), 637-669
%-------------------------------------------------------------

endogenous C Z K

exogenous E

parameters alpha delta beta upsilon s_tp_1_2 s_tp_2_1

parameters(s,2) mu sigma rho

model

	C^(upsilon-1) = beta*Z^(upsilon-1)*C{+1}^(upsilon-1)*(alpha*Z{+1}^(1-alpha)*K^(alpha-1)+1-delta);

	C + Z*K = Z^(1-alpha)*K{-1}^alpha+(1-delta)*K{-1};

	log(Z) = (1-rho)*mu +rho*log(Z{-1})+sigma*E;

steady_state_model

	Z=exp(mu);
	
	K=(1/alpha*exp((alpha-1)*mu)*(1/beta*exp((1-upsilon)*mu)-1+delta))^(1/(alpha-1));
	
	C = K*(1-delta-exp(mu)+1/alpha*(1/beta*exp((1-upsilon)*mu)-1+delta));

parameterization
	alpha, 0.33;
	beta, 0.9976;
	upsilon,-1;
	delta, 0.025;
	mu(s,1),0.0274;
	mu(s,2),-0.0337;
	rho(s,1),0.1;
	rho(s,2),0.0;
	sigma(s,1), 0.0072;
	sigma(s,2), 0.0216;
	s_tp_1_2, 1-0.75;
	s_tp_2_1, 1-0.5;
	 