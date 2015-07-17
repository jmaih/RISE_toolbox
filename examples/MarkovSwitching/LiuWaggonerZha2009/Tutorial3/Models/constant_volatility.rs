parameters sig_pai, "$\sigma_{\pi}$",	sig_y, "$\sigma_{y}$", sig_i, "$\sigma_{i}$"

parameterization
	sig_pai,   		0.1, 	0.05, 	  3, 	weibull_pdf(0.9); 
	sig_y,  		0.1, 	0.05, 	  3, 	weibull_pdf(0.9);  
	sig_i, 			0.1, 	0.05,     3, 	weibull_pdf(0.9); 

	