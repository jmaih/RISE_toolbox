% In passing, we demonstrate a little bit of the macro language in RISE
% we include the previous file that does not contain the steady_state_model block

@# include "fs2000.rs"


steady_state_model
	@# if approximate==1
		k = 6;
		m = mst;
		P = 2.25;
		c = 0.45;
		e = 1;
		W = 4;
		R = 1.02;
		d = 0.85;
		n = 0.19;
		l = 0.86;
		y = 0.6;
		gy_obs = gam;  %	gy_obs = exp(gam);
		gp_obs = -gam; %	gp_obs = exp(-gam);
		dA = exp(gam);
	@# else
	%NB: xx_ssmdef_1 ... xx_ssmdef_9 are known words to rise. They are used to define auxiliary
	% constants in the steady_state_model block
		dA = exp(gam);
		xx_ssmdef_1 = 1/dA;
		m = mst;
		xx_ssmdef_2 = ( (1-xx_ssmdef_1*bet*(1-del)) / (alp*xx_ssmdef_1^alp*bet) )^(1/(alp-1));
		xx_ssmdef_3 = ( ((xx_ssmdef_2*xx_ssmdef_1)^alp - (1-xx_ssmdef_1*(1-del))*xx_ssmdef_2)/mst )^(-1);
		xx_ssmdef_4 = psi*mst^2/( (1-alp)*(1-psi)*bet*xx_ssmdef_1^alp*xx_ssmdef_2^alp );
		n  = xx_ssmdef_3/(xx_ssmdef_4+xx_ssmdef_3);
		P  = xx_ssmdef_3 + xx_ssmdef_4;
		k  = xx_ssmdef_2*n;
		l  = psi*mst*n/( (1-psi)*(1-n) );
		c  = mst/P;
		d  = l - mst + 1;
		y  = k^alp*n^(1-alp)*xx_ssmdef_1^alp;
		R  = mst/bet;
		W  = l/n;
		e = 1;
		gp_obs = log(m/dA);% accommodating dynare's loglinear option
		gy_obs = log(dA); % accommodating dynare's loglinear option
	@# end
