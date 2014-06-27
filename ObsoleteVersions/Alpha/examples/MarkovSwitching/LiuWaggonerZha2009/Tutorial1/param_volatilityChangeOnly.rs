% in this case, only sigd is controlled by markov chain nk
% the other parameters are controlled by the const markov chain

parameters tau rhod rhos sigs sigr kappa beta rhor gamma_1 gamma_2

parameters(nk,2) sigd

parameterization
	tau          ,   0.6137;
	rhod         ,   0.7550;
	rhos         ,   0.835 ;
	sigd(nk,1)   ,   0.2250;
	sigd(nk,2)   ,   5*0.2250;
	sigs         ,   0.3712;
	sigr         ,   0.2050;
	kappa        ,   0.6750;
	beta         ,   0.9949;
	rhor         ,   0.72  ;
	gamma_1      ,   2.19  ;
	gamma_2      ,   0.235 ;
	