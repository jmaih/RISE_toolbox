% financial frictions (page 34)
endogenous	X R P G U I E

exogenous EG EU	EE

parameters beta kappa sigu sigg sige rhou rhog rhoe chi	psi

observables P I R

model

	P = beta*P{+1}+kappa*X;
	
	X = X{+1}-(I-P{+1}-G);

	R = psi*P+U;

	I = chi*R+E;

	G = rhog*G{-1} + sigg*EG;

	U = rhou*U{-1} + sigu*EU;

	E = rhoe*E{-1} + sige*EE;

	

