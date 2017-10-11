% Including a second lag of a control variable (page 56)
endogenous C N W R

exogenous  ER EW

parameters h gam rhow sigw sigr b1

observables N C

model

	N = b1*N{-2} + W - gam*C;

	C = (1-h)*W+h*C{+1}+R;

	W = rhow*W{-1}+sigw*EW;
	
	R = sigr*ER;
	