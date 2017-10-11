% Including an observed exogenous variable (page 59)
endogenous	X R P G U E

exogenous EG EU	EE

parameters beta kappa sigu sigg rhou rhog sige rhoe	psi

observables R P E

model

	X = X{+1}-(R-P{+1}-G);

	R = 1/beta*P+U;

	P = beta*P{+1}+kappa*X+psi*E;
	
	U = rhou*U{-1} + sigu*EU;

	G = rhog*G{-1} + sigg*EG;

	E = rhoe*E{-1} + sige*EE;
