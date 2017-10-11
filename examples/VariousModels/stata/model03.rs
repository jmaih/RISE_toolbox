% New classical model (page 27-28)
endogenous	C R H W X K Z G	Y

exogenous EG EZ

parameters beta delta phi1 phi2 eta alpha rhoz rhog sigz sigg

observables  Y 

model

	C = C{+1}-(1-beta+beta*delta)*R{+1};

	eta*H = W - C;
	
	phi1*X = Y -phi2*C -G;

	Y = (1-alpha)*(Z+H)+alpha*K{-1};

	W = Y - H;

	R = Y-K{-1};

	K = delta*X + (1-delta)*K{-1};

	Z = rhoz*Z{-1} + sigz*EZ;

	G = rhog*G{-1} + sigg*EG;

	

