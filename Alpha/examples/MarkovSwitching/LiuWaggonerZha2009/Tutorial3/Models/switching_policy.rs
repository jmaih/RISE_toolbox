parameters pol_tp_1_2 pol_tp_2_1

parameters(pol,2) rho_i, "$\rho_{i}$" gam_y, "$\gamma_{y}$" gam_pai, "$\gamma_{\pi}$" c_i, "$c_{i}$"

parameterization
	rho_i(pol,1), 			0.6, 	0.1, 	0.7, 	beta_pdf(0.9);
	rho_i(pol,2), 			0.6, 	0.1, 	0.7, 	beta_pdf(0.9);
	gam_y(pol,1), 			0.5, 	0.1, 	1.5, 	gamma_pdf(0.9);
	gam_y(pol,2), 			0.5, 	0.1, 	1.5, 	gamma_pdf(0.9);
	gam_pai(pol,1), 		1.5, 	0.5, 	3, 		gamma_pdf(0.9);
	gam_pai(pol,2), 		1.0, 	0.5, 	3, 		gamma_pdf(0.9);
	c_i(pol,1), 			0  , 	-1, 	1, 		normal_pdf(0.9); 
	c_i(pol,2), 			0  , 	-1, 	1, 		normal_pdf(0.9);
	% transition probabilities
	pol_tp_1_2,   			0.15,	0.1, 	0.5, 	beta_pdf(0.9);
	pol_tp_2_1,   			0.15,	0.1, 	0.5, 	beta_pdf(0.9);

parameter_restrictions
	% high monetary policy response in the first regime
	gam_pai(pol,1)>=gam_pai(pol,2);	
