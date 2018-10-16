%-------------------------------------------------------------
% Nonlinear New Keynesian Model with habits
% Reference: Foerster, Rubio-Ramirez, Waggoner and Zha (2013)
% Perturbation Methods for Markov Switching Models.
%-------------------------------------------------------------
% TODO: reinstate the possibility of having leads and lags on parameters?

endogenous PAI, "Inflation", C, "Consumption", R, "Interest rate", LAMBDA "Lagrange Multiplier", MU "Auxiliary variable"

exogenous E	"Monetary policy shock"

parameters beta, kappa, eta, varphi, sigma, a_tp_1_2, a_tp_2_1

parameters(a,2) mu, psi

model

	LAMBDA = 1/(C-varphi*exp(-mu)*C{-1})-beta*varphi/(C{+1}*exp(MU{+1})-varphi*C);

	LAMBDA = beta*LAMBDA{+1}/exp(MU{+1})*R/PAI{+1};

	kappa*(PAI-1)*PAI = (1-eta)+eta/LAMBDA+beta*kappa*(PAI{+1}-1)*PAI{+1}*LAMBDA{+1}/LAMBDA*C{+1}/C*(1-kappa/2*(PAI-1)^2)/(1-kappa/2*(PAI{+1}-1)^2);
	
	R = R{stst}*PAI^psi*exp(sigma*E);

	MU=mu;

steady_state_model
    R=exp(mu)/beta;
    PAI=1;
	LAMBDA=eta/(eta-1);
    C=(exp(mu)-beta*varphi)/(exp(mu)-varphi)*(eta-1)/eta;
	MU=mu;

parameterization
	beta, .9976;
	kappa, 161;
	eta, 10;
	varphi, 0.7;
	sigma, 0.0025;
	mu(a,1), 0.005+0.0025;
	mu(a,2), 0.005-0.0025;
	psi(a,1), 3.1;
	psi(a,2), 0.9;
	a_tp_1_2,1-.9; 
	a_tp_2_1,1-.9;
