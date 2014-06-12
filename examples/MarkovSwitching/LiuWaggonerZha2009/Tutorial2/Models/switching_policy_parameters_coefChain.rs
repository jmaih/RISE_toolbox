% We declare the parameters of the markov chain that will control the
% coefficients. The parameters of that markov chain are themselves
% controlled by the const markov chain since they are constant.
parameters coef_tp_1_2 coef_tp_2_1

parameters(coef,2) gamma_1, "$\gamma _{1}$", gamma_2, "$\gamma _{2}$"

parameterization
	gamma_1(coef,1) ,    2.19  ,     0.5000,    5.0000,  gamma_pdf(.90);
	gamma_1(coef,2) ,    0.77  ,     0.5000,    5.0000,  gamma_pdf(.90);
	gamma_2(coef,1) ,    0.30  ,     0.0500,    3.0000,  gamma_pdf(.90);
	gamma_2(coef,2) ,    0.17  ,     0.0500,    3.0000,  gamma_pdf(.90);
	% transition probabilities
	coef_tp_1_2     ,   0.0128 ,    0.0500,    0.1500,   beta_pdf(.90);
	coef_tp_2_1     ,   0.0128 ,    0.0500,    0.1500,   beta_pdf(.90);

% for identification purposes we need to impose the regime in which a particular
% parameter will assume the greatest value. We choose to identify the first
% regime as the regime with the highest response of the

parameter_restrictions
	gamma_1(coef,1)>=gamma_1(coef,2)

