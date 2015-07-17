% we declare the endogenous variables
endogenous PAI, R

% we declare the exogenous variables
exogenous E
% All shocks have standard deviation 1. In this sense
% there is no such a thing as a covariance matrix
% this is by pure convenience

% we declare the constant parameters
parameters a_tp_1_2, a_tp_2_1 
% we declare the switching parameters
parameters(a,2) betta, delta, phi, rho

% This is the model block. The linear statement in parenthesis is useless
% the parser won't even read it. I just leave it there to remind people that
% the procedure will not do linearization for you as more elaborate procedures
% like dynare would.
model(linear)
	phi*PAI=PAI(+1)+delta*PAI(-1)+betta*R;
	R=rho*R(-1)+E;