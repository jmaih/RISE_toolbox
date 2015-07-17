% "The inflation bias under Calvo and Rotemberg pricing",
% by Campbell Leith and Ding Liu. Working paper 2014 06,
% University of Glasgow.

endogenous R, DELTA, Z, C, W, N, K, J, PSTAR, PI, MC

exogenous EZ, ER

parameters Beta, sigma, psi, epsilon, theta, kappa_c, kappa_pi, rho_r, rho, phi,
std_EZ	
std_ER

model

    Beta*R*((C/C(1))^sigma*(1/PI(1))) - 1;
	
    W - phi*N^psi*C^sigma;
	
    C - Z*N/DELTA;
	
    PSTAR - (epsilon/(epsilon - 1))*(K/J);
	
    K - C^(1-sigma)*MC - theta*Beta*(PI(1))^epsilon*K(1);
	
    J - C^(1-sigma) - theta*Beta*PI(1)^(epsilon-1)*J(1);
	
    1 - (1 - theta)*PSTAR^(1-epsilon) - theta*PI^(epsilon-1);
	
    DELTA - (1-theta)*PSTAR^(-epsilon) - theta*(PI)^epsilon*DELTA(-1);
	
    log(R/steady_state(R)) - (1-rho_r)*(kappa_c*log(C/steady_state(C)) +
		kappa_pi*log(PI/steady_state(PI))) - rho_r*log(R(-1)/steady_state(R)) - std_ER*ER;
	
    MC - W*N/C;
	
    log(Z/steady_state(Z)) - rho*log(Z(-1)/steady_state(Z)) - std_EZ*EZ;

steady_state_model
	PI = 1;
	DELTA = 1;
	Z = 1;
	N = 1/3;
	C = Z*N;
	MC = (epsilon-1)/epsilon;
	W = C*MC/N;
	phi = W/(N^psi*C^sigma);
	K = C^(1-sigma)*MC/(1 - theta*Beta*PI^epsilon);
	J = C^(1-sigma)/(1 - theta*Beta*PI^(epsilon-1));
	PSTAR = (epsilon/(epsilon - 1))*(K/J);
	R = PI/Beta;

parameterization
	Beta , 0.99;
	sigma , 2;
	psi , 3;
	epsilon , 6;
	theta , 0.75;
	kappa_c , 0.5;
	kappa_pi , 1.5;
	rho_r , 0.8;
	rho , 0.8;
	std_EZ, 0.01;
	std_ER, 0.01;
