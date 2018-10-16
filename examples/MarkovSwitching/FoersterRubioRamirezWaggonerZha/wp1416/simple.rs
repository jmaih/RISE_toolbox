%-------------------------------------------------------------
% Simple inflation model
% Reference: Foerster, Rubio-Ramirez, Waggoner and Zha (2014)
% Perturbation Methods for Markov Switching Models.
% Atlanta Fed working paper 2014-16
%-------------------------------------------------------------

endogenous PAI

exogenous E

parameters s_tp_1_2 s_tp_2_1

parameters(s,2) phi sigma

model

	phi*PAI+sigma*E=PAI{+1};

parameterization
    phi(s,1),1.25;
	phi(s,2),0.96;
	sigma(s,1), 0.1;
	sigma(s,2), 0.6;
	s_tp_1_2, 1-0.95;
	s_tp_2_1, 1-0.85;
	 