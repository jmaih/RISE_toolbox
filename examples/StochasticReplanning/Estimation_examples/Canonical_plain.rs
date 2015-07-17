endogenous	DPQ_P_NW, Y, ZGDP, ZI, ZPAI, ZY, D_GDP_NW, I, PAI, R, RN3M_NW

exogenous EGDP, EI, EPAI, EY

parameters beta_lag, beta_lead, beta_r, gam_lag, gam_y, gyss, iss, lamb_lag,
lamb_lead, lamb_y, paiss, rhogdp, rhoi, rhopai, rhoy, siggdp, sigi, sigpai, sigy

parameters prob_commit

model(linear)
   
   Y=beta_lag*Y(-1)+beta_lead*Y(+1)-beta_r*R(-1)+ZY;

   PAI=lamb_lag*PAI(-1)+lamb_lead*PAI(+1)+lamb_y*Y(-1)+ZPAI;

//   // Taylor rule
//   I=gam_lag*I(-1)+(1-gam_lag)*(PAI(+4)+gam_y*Y)+ZI;

   R=I-PAI(+1);

   D_GDP_NW=Y-Y(-1)+ZGDP;

   DPQ_P_NW=paiss+PAI;

   RN3M_NW=iss+I;

   ZI=rhoi*ZI(-1)+sigi*EI;
   
   ZPAI=rhopai*ZPAI(-1)+sigpai*EPAI;
   
   ZY=rhoy*ZY(-1)+sigy*EY;
   
   ZGDP=(1-rhogdp)*gyss+rhogdp*ZGDP(-1)+siggdp*EGDP;

// the definition in the model is used in specifying the planner objective
planner_objective{discount = 0.99,commitment=prob_commit} -.5*(1*PAI^2+.6*Y^2);	 

observables DPQ_P_NW, D_GDP_NW, RN3M_NW;

parameterization
// not estimated
	gyss   		 ,0 	;						 
	iss    		 ,0 	;						 
	paiss  		 ,0 	;
	beta_lag 	 ,0.5000;     
	beta_lead	 ,0.4000;	 
	beta_r  	 ,0.9000;	 
	gam_lag 	 ,0.6000;	 
	gam_y   	 ,0.5000;	 
	lamb_lag	 ,0.8000;	 
	lamb_lead    ,0.1000;	 
	lamb_y  	 ,0.3000;	 
	rhogdp 		 ,0.5000;	 
	rhoi   		 ,0.5000;	 
	rhopai 		 ,0.5000;	 
	rhoy   		 ,0.5000;     
	siggdp 		 ,0.0050, 0.0001, 0.5000;	 		 
	sigi   		 ,0.0050, 0.0001, 0.5000;	 
	sigpai 		 ,0.0050, 0.0001, 0.5000;	 
	sigy   		 ,0.0050, 0.0001, 0.5000;
	// probability of commitment 
	prob_commit     ,0.5000, 0.0000, 1.0000;
