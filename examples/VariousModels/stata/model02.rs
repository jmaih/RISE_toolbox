% New Keynesian model (page 22)
endogenous	X R P G U

exogenous EG EU

parameters beta kappa sigu sigg rhou rhog psi

observables  P R

model

	P = beta*P{+1}+kappa*X;
	
	X = X{+1}-(R-P{+1}-G);

	R = psi*P+U;

	U = rhou*U{-1} + sigu*EU;

	G = rhog*G{-1} + sigg*EG;

	

