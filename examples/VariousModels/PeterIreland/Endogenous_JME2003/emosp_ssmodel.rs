@#include "emosp.rs"

steady_state_model
	X=1;
	A=1;
	E=e_ss;
	Z=z_ss;
	V=1;
	% mu determined by policy
	MU=mu_ss;
	PAI=MU/g;
	R=MU/beta;
	Q=g/beta-1+delta;
	% lambda comes here
	LAMBDA=(
		eta+(1-alpha)*1/(1+E*(R/(R-1))^(gam-1))*1/(theta/(theta-1)-(g-1+delta)*alpha/Q)
		)
		/(
		(1-alpha)*Z*((theta-1)/theta)^(1/(1-alpha))*(alpha/Q)^(alpha/(1-alpha))
		);
	XI = (theta-1)/theta*LAMBDA;
	C = 1/(1+E*(R/(R-1))^(gam-1))*(1/LAMBDA);
	M = E*(R/(R-1))^gam*C;
	Y = 1/(1-(g-1+delta)*(theta-1)/theta*alpha/Q)*C;
	K = (theta-1)/theta*alpha*Y/Q;
	I = (g-1+delta)*K;
	H=1/Z*(Y/K^alpha)^(1/(1-alpha));
	W=(1-alpha)*(theta-1)/theta*Y/H;
	D = Y-W*H-Q*K;
	N=Y/H;
	LC=log(C);
	LI=log(I);
	LM=log(M);
	LPI=log(PAI);
	LR=log(R);	