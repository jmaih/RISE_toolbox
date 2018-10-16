%-------------------------------------------------------------
% Nonlinear Nonstationary RBC model
% Reference: Foerster, Rubio-Ramirez, Waggoner and Zha (2013)
% Perturbation Methods for Markov Switching Models.
%-------------------------------------------------------------

endogenous C Z K

exogenous E

parameters alpha delta sigma beta s_tp_1_2 s_tp_2_1

parameters(s,2) mu

model

	1/C = beta*1/C{+1}*(alpha*Z{+1}*K^(alpha-1)+1-delta);

	C + K = Z*K{-1}^alpha+(1-delta)*K{-1};

	log(Z) = mu+log(Z{-1})+sigma*E;
	 