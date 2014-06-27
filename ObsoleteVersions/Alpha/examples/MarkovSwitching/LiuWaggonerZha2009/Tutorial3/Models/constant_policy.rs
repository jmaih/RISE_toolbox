parameters rho_i, "$\rho_{i}$", gam_y, "$\gamma_{y}$", gam_pai, "$\gamma_{\pi}$", c_i, "$c_{i}$"

parameterization
	rho_i, 			0.6, 	0.1, 	0.7, 	beta_pdf(0.9);
	gam_y, 			0.5, 	0.1, 	1.5, 	gamma_pdf(0.9);
	gam_pai, 		1.5, 	0.5, 	3, 		gamma_pdf(0.9);
	c_i, 			0  , 	-1, 	1, 		normal_pdf(0.9); 
