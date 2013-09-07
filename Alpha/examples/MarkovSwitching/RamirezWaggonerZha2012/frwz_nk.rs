%-------------------------------------------------------------
% Nonlinear New Keynesian Model
% Reference: Foerster, Rubio-Ramirez, Waggoner and Zha (2012)
% Perturbation Methods for Markov Switching Models.
%-------------------------------------------------------------
endogenous PAI,Y,R

exogenous EPS_R

parameters a_tp_1_2, a_tp_2_1, betta, eta, kappa, mu_bar, rhor sigr

parameters(a,2) mu, psi

model
	1-betta*(1-.5*kappa*(PAI-1)^2)*Y*R/((1-.5*kappa*(PAI{+1}-1)^2)*Y{+1}*exp(mu)*PAI{+1});
	
	1-eta+eta*(1-.5*kappa*(PAI-1)^2)*Y+betta*kappa*(1-.5*kappa*(PAI-1)^2)*(PAI{+1}-1)*PAI{+1}/(1-.5*kappa*(PAI{+1}-1)^2)
	-kappa*(PAI-1)*PAI;

	(R{-1}/steady_state(R))^rhor*(PAI/steady_state(PAI))^((1-rhor)*psi)*exp(sigr*EPS_R)-R/steady_state(R);

steady_state_model(unique)%imposed
    PAI=1;
    Y=(eta-1)/eta;
    R=exp(mu_bar)/betta*PAI;

parameterization
	a_tp_1_2,1-.9; 
	a_tp_2_1,1-.9;
	betta, .99;
	kappa, 161;
	eta, 10;
	rhor, .8;
	sigr, 0.0025;
	mu_bar,0.02; 
	mu(a,1), 0.03;
	mu(a,2), 0.01;
	psi(a,1), 3.1;
	psi(a,2), 0.9;
