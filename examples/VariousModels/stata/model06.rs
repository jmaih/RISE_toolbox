%Including a lag on a control variable (page 47)
endogenous P Y R Z U

exogenous EZ EU

parameters beta kappa rhor rhoz sigz sigu rhou

observables P R

model

	P = beta*P{+1}+kappa*Y;

	Y = Y{+1} - (R-P{+1}-rhoz*Z);

	R = rhor*R{-1}+(1-rhor)*(1/beta*P+U);

	Z = rhoz*Z{-1}+sigz*EZ;

	U = rhou*U{-1}+sigu*EU;