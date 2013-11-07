% We declare the parameters of the markov chain that will control the
% volatility of the shocks. The parameters of that markov chain are themselves
% controlled by the const markov chain since they are constant.
parameters vol_tp_1_2, vol_tp_2_1

parameters(vol,2) sigd, "$\sigma _{d}$" sigs, "$\sigma _{s}$"	sigr, "$\sigma _{r}$"

parameterization
	sigd(vol,1)          ,    0.18   ,    0.0005,    1.0000,  weibull_pdf(.90);
	sigd(vol,2)          ,    0.27   ,    0.0005,    1.0000,  weibull_pdf(.90);
	sigs(vol,1)          ,    0.3712 ,    0.0005,    1.0000,  weibull_pdf(.90);
	sigs(vol,2)          ,    0.8701 ,    0.0005,    1.0000,  weibull_pdf(.90);
	sigr(vol,1)          ,    0.18   ,    0.0005,    1.0000,  weibull_pdf(.90);
	sigr(vol,2)          ,    0.23   ,    0.0005,    1.0000,  weibull_pdf(.90);
	% transition probabilities
	vol_tp_1_2           ,   0.0128  ,    0.0500,    0.1500,  beta_pdf(.90);
	vol_tp_2_1           ,   0.0128  ,    0.0500,    0.1500,  beta_pdf(.90);

% for identification purposes we need to impose the regime in which a particular
% parameter will assume the greatest value. We choose to identify the second
% regime as the regime with the highest volatility

parameter_restrictions
	sigd(vol,2)>=sigd(vol,1);
