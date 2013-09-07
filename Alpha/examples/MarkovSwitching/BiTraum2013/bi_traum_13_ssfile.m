function bi_traum_13_ssfile()

	shat = stilde - 0.6*steady_state(Y);

	eta2 = 1/(stilde-shat)*log(ptilde/phat*(1-phat)/(1-ptilde));

	eta1 = log(ptilde/(1-ptilde)) - eta2*stilde;


end