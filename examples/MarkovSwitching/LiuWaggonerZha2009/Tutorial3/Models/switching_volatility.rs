parameters vol_tp_1_2, vol_tp_2_1

parameters(vol,2) sig_pai, "$\sigma_{\pi}$" sig_y, "$\sigma_{y}$", sig_i, "$\sigma_{i}$" 

parameterization
	sig_pai(vol,1),   		0.1, 	0.05, 	  3, 	weibull_pdf(0.9); 
	sig_pai(vol,2),   		0.1, 	0.05, 	  3, 	weibull_pdf(0.9); 
	sig_y(vol,1),  		    0.1, 	0.05, 	  3, 	weibull_pdf(0.9);  
	sig_y(vol,2),  		    0.1, 	0.05, 	  3, 	weibull_pdf(0.9);  
	sig_i(vol,1), 			0.1, 	0.05,     3, 	weibull_pdf(0.9); 
	sig_i(vol,2), 			0.1, 	0.05,     3, 	weibull_pdf(0.9); 
	% transition probabilities
	vol_tp_1_2,   			0.15,	0.1, 	0.5, 	beta_pdf(0.9);
	vol_tp_2_1,   			0.15,	0.1, 	0.5, 	beta_pdf(0.9);

parameter_restrictions
	% high volatility is the second regime
	sig_pai(vol,2)>=sig_pai(vol,1);
