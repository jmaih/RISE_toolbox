endogenous P Y R U Z

exogenous  EZ EU

parameters beta gam sigu sigz rhoz kappa rhou

observables P R

model

	P = beta*P{+1}+kappa*Y;

	Y = Y{+1}-gam*(R-P{+1}-rhoz*Z);

	beta*R = P+beta*U;

	Z=rhoz*Z{-1}+sigz*EZ;

	U=rhou*U{-1}+sigu*EU;