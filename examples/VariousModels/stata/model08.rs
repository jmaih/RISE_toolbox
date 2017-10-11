% Including an expectation of a control dated by more than one period (page 53)
endogenous C N W R

exogenous  ER EW

parameters h gam rhow rhor sigw sigr

observables C N

model

	C = (1-h)*W+h*C{+2}+R;

	N = W - gam*C;

	W = rhow*W{-1}+sigw*EW;
	
	R = rhor*R{-1}+sigr*ER;
	