% Specifying a shock on a control variable (page 45)
endogenous C N W

exogenous  E EW

parameters h sige gam rho sigw

observables C N

model

	C = (1-h)*W+h*C{+1}+sige*E;

	N = W - gam*C;

	W = rho*W{-1}+sigw*EW;
	