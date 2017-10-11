% Correlated state variables (page 62)
endogenous	Y P Z G

exogenous EG EZ

parameters alpha rhog rhogz rhoz sigg sigz

observables Y P

model

	Y = Y{+1}+alpha*P+G;

	P = Z;

	G = rhog*G{-1} + rhogz*Z{-1}+sigg*EG;

	Z = rhoz*Z{-1} + sigz*EZ;
