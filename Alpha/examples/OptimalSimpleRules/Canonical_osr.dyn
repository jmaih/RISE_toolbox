var	Y, I, PAI, R, DI;

varexo EI, EPAI, EY;

parameters beta_lag, beta_lead, beta_r, gam_lag, gam_y, lamb_lag,
lamb_lead, lamb_y, sigi, sigpai, sigy;

model(linear);

   Y=beta_lag*Y(-1)+beta_lead*Y(+1)-beta_r*R(-1)+sigy*EY;

   PAI=lamb_lag*PAI(-1)+lamb_lead*PAI(+1)+lamb_y*Y(-1)+sigpai*EPAI;

   // Taylor rule
   I=gam_lag*I(-1)+(1-gam_lag)*(PAI(+4)+gam_y*Y)+sigi*EI;
				  
   R=I-PAI(+1);
   // interest rate smoothing
   DI=I-I(-1);
end;

// no need to specify a list of observable variables

// specify the planner's objective. 
planner_objective -.5*(1*PAI^2+.3*Y^2+0.9*DI^2);//{discount = 0.99}
// The discount factor does not play a role. So I can... neglect it.
//planner_objective{discount = 0.99} -.5*(1*PAI^2+.3*Y^2+0.9*DI^2);	 

parameterization;
// not estimated
	beta_lag 	 ,0.5000;     
	beta_lead	 ,0.4000;	 
	beta_r  	 ,0.9000;
	// estimate the monetary policy parameter
	gam_lag 	 ,0.6000, 0.0000, 1.0000;	 
	gam_y   	 ,0.5000, 0.0000, 10.0000;	 
	lamb_lag	 ,0.8000;	 
	lamb_lead    ,0.1000;	 
	lamb_y  	 ,0.3000;	 
	sigi   		 ,0.5000;	 
	sigpai 		 ,0.5000;	 
	sigy   		 ,0.5000;	 
end;
