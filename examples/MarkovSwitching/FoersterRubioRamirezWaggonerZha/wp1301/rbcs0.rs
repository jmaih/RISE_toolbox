%-------------------------------------------------------------
% Nonlinear stationary RBC model
% Reference: Foerster, Rubio-Ramirez, Waggoner and Zha (2013)
% Perturbation Methods for Markov Switching Models.
%-------------------------------------------------------------

endogenous C Z K

exogenous E

parameters alpha delta sigma beta mu

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
	mu,0.00333+0.00167;	 