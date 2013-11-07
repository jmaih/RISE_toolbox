parameters(vol,2) gamma_1, "$\gamma _{1}$", gamma_2, "$\gamma _{2}$"

parameterization
	gamma_1(vol,1) ,    2.19  ,     0.5000,    5.0000,  gamma_pdf(.90);
	gamma_1(vol,2) ,    0.77  ,     0.5000,    5.0000,  gamma_pdf(.90);
	gamma_2(vol,1) ,    0.30  ,     0.0500,    3.0000,  gamma_pdf(.90);
	gamma_2(vol,2) ,    0.17  ,     0.0500,    3.0000,  gamma_pdf(.90);

