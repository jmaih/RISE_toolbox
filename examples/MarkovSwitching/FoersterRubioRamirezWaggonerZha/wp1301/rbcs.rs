%-------------------------------------------------------------
% Nonlinear stationary RBC model
% Reference: Foerster, Rubio-Ramirez, Waggoner and Zha (2013)
% Perturbation Methods for Markov Switching Models.
%-------------------------------------------------------------

endogenous C Z K

exogenous E

parameters alpha delta sigma beta s_tp_1_2 s_tp_2_1

parameters(s,2) mu

model

	1/C = beta*Z^(1/(alpha-1))/C{+1}*(alpha*Z{+1}*K^(alpha-1)+1-delta);

	C + K*Z^(1/(alpha-1)) = Z*K{-1}^alpha+(1-delta)*K{-1};

	log(Z) = mu + sigma*E;

steady_state_model

	Z=1;
	
	K=(1/(alpha*exp(mu))*(1/(beta*exp(mu/(alpha-1)))-1+delta))^(1/(alpha-1));
	
	C = exp(mu)*K^alpha+(1-delta)*K-K*exp(mu/(1-alpha));

parameterization
	alpha, 0.33;
	beta, 0.9976;
	delta, 0.025;
	sigma, 0.0002;
	mu(s,1),0.00333+0.00167;
	mu(s,2),0.00333-0.00163;
	s_tp_1_2, 1-0.9;
	s_tp_2_1, 1-0.9;
	 