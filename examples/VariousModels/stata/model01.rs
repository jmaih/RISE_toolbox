% Specifying a DSGE model (page 6-7)
endogenous	X R P G U

exogenous EG EU

parameters beta kappa sigu sigg rhou rhog

observables P R

model

	X = X{+1}-(R-P{+1}-G);

	P = beta*P{+1}+kappa*X;
	
	R = 1/beta*P+U;

	U = rhou*U{-1} + sigu*EU;

	G = rhog*G{-1} + sigg*EG;

	

