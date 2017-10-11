%Including a lag of a state variable (page 50)
endogenous P Y R Z U

exogenous EZ EU

parameters beta kappa rhoz1 rhoz2 sigz sigu rhou

observables P R

model

	P = beta*P{+1}+kappa*Y;

	Y = Y{+1} - (R-P{+1}-Z);

	R = 1/beta*P+U;

	Z = rhoz1*Z{-1}+rhoz2*Z{-2}+sigz*EZ;

	U = rhou*U{-1}+sigu*EU;