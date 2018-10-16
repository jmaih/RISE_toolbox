%-------------------------------------------------------------
% Nonlinear New Keynesian Model
% Reference: Foerster, Rubio-Ramirez, Waggoner and Zha (2013)
% Perturbation Methods for Markov Switching Models.
%-------------------------------------------------------------

endogenous PAI, "Inflation", Y, "Output gap", R, "Interest rate"

exogenous EPS_R	"Monetary policy shock"

parameters betta, eta, kappa, rhor sigr a_tp_1_2, a_tp_2_1

parameters(a,2) mu, psi

model
	1=betta*(1-.5*kappa*(PAI-1)^2)*Y*R/((1-.5*kappa*(PAI{+1}-1)^2)*Y{+1}*exp(mu)*PAI{+1});
	
	1-eta+eta*(1-.5*kappa*(PAI-1)^2)*Y+betta*kappa*(1-.5*kappa*(PAI-1)^2)*(PAI{+1}-1)*PAI{+1}/(1-.5*kappa*(PAI{+1}-1)^2)
	-kappa*(PAI-1)*PAI;

	(R{-1}/R{stst})^rhor*(PAI/PAI{stst})^((1-rhor)*psi)*exp(sigr*EPS_R)-R/R{stst};

steady_state_model
    PAI=1;
    Y=(eta-1)/eta;
    R=exp(mu)/betta*PAI;

parameterization
	a_tp_1_2,1-.9; 
	a_tp_2_1,1-.9;
	betta, .9976;
	kappa, 161;
	eta, 10;
	rhor, .8;
	sigr, 0.0025;
	mu(a,1), 0.005+0.0025;
	mu(a,2), 0.005-0.0025;
	psi(a,1), 3.1;
	psi(a,2), 0.9;
