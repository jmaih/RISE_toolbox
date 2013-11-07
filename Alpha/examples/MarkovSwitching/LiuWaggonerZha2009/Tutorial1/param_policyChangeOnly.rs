% in this case only gamma_1 changes. And so,
% gamma_1 is controlled by markov chain nk
% while the other parameters are controlled by
% the constant markov chain
parameters tau rhod rhos sigd sigs sigr kappa beta rhor gamma_2

parameters(nk,2) gamma_1

parameterization
	tau          ,   0.6137;
	rhod         ,   0.7550;
	rhos         ,   0.835 ;
	sigd         ,   0.2250;
	sigs         ,   0.6206;
	sigr         ,   0.2050;
	kappa        ,   0.6750;
	beta         ,   0.9949;
	rhor         ,   0.72  ;
	gamma_1(nk,1),   2.19  ;
	gamma_1(nk,2),   0.77  ;
	gamma_2      ,   0.235 ;
