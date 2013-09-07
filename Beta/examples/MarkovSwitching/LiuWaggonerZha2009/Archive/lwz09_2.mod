// Zheng Liu, Daniel F. Waggoner and Tao Zha (2009):
//Asymmetric expectation effects of regime shifts in monetary policy
// Review of Economic Dynamics 12 (2009) 284-303

// this model: page 290
// parameterization not given in the paper...

var A PAI $Inflation$ V MU R $Interest rate$ Y $Output$;

varexo
EPS_A, $Preference$	// preference shock
EPS_MU, $Markup$	// Markup shock
EPS_R, $Monetary Policy$   // Monetary policy shock
EPS_V, $Technology$;	// technology shock

parameters alpha,  b,  beta, eta, iota, lambda, phipai, phiy, q_tp_1_2, q_tp_2_1, rhoa, rhomu, rhor, rhov, siga
sigmu, sigr, sigv, theta, xi;

model
	# psi=(1-beta*eta)*(1-eta)/eta*1/(1+theta*(1-alpha)/alpha);	// definition page 289
	
	A=rhoa*A{-1}+siga*EPS_A; // eq 1
	
	PAI-iota*PAI{-1}=beta*(PAI{+1}-iota*PAI)+psi*((xi+1)/alpha*Y+b/(lambda-b)*(Y-Y{-1}+V))+psi*MU; // eq 16

	V=rhov*V{-1}+sigv*EPS_V; // eq 17

	MU=rhomu*MU{-1}+sigmu*EPS_MU; // eq 18

	Y{+1}-(lambda+b)/lambda*Y+b/lambda*Y{-1} = (1-b/lambda)*(R-PAI{+1})+(b/lambda-rhov)*V
		-(lambda-b)*(1-rhoa)/lambda*A; // eq 19
	
	R=rhor*R{-1}+(1-rhor)*(phipai*PAI+phiy*Y)+sigr*EPS_R; // eq 20

end

parameterization
	// Preferences
	beta       , 0.9952;
	xi         , 2.0000;
	b          , 0.7500; // habit formation
	// Technology
	alpha      , 0.7000; // elasticity of output wrt labor
	lambda     , 1.0050;
	theta      , 10.000; // elasticity of subst. btw interm goods
	// Price setting
	eta(q,1)   , 0.6600; // Calvo stickiness
	eta(q,2)   , 0.6600; 
	iota(q,1)  , 1.0000; // inflation indexation
	iota(q,2)  , 1.0000; 
	// policy rule
	rhor       , 0.5500;
	phiy       , 0.5000;
	// Aggregate shocks
	rhoa       , 0.9000;
	rhomu      , 0.9000;
	rhov       , 0.0000;
	siga       , 0.1000;
	sigmu      , 0.1000;
	sigr       , 0.1000;
	sigv       , 0.1000;
	// Regime transition probability
	q_tp_1_2   , 1-0.95;
	q_tp_2_1   , 1-0.95;
	phipai(q,1), 0.9000; // dovish
	phipai(q,2), 2.5000; // hawkish
end