endogenous	DPQ_P_NW "Inflation", Y, ZGDP, ZI, ZPAI, ZY, D_GDP_NW "Growth", I, PAI, R, RN3M_NW "Fed Funds Rate"

exogenous EGDP "output shock",EI "monetary policy shock",EPAI "Cost push shock",EY "IS shock"

parameters beta_lag "\beta_{lag}", beta_lead "\beta_{lead}", beta_r "\beta_{r}",
gam_lag "\gamma_{lag}", gam_y "\gamma_{y}", gyss, iss, lamb_lag "\lambda_{lag}",
lamb_lead "\lambda_{lead}", lamb_y "\lambda_{y}", paiss,
rhogdp "\rho_{gdp}", rhoi "\rho_{i}", rhopai "\rho_{\pi}", rhoy "\rho_{y}",
siggdp "\sigma_{gdp}", sigi "\sigma_{i}", sigpai "\sigma_{\pi}", sigy "\sigma_{y}"

model

	# junk=beta_lead;
	
   Y=beta_lag*Y(-1)+beta_lead*Y(+1)-beta_r*R(-1)+ZY;

   PAI=lamb_lag*PAI(-1)+lamb_lead*PAI(+1)+lamb_y*Y(-1)+ZPAI;

   I=gam_lag*I(-1)+(1-gam_lag)*(PAI(+4)+gam_y*Y)+ZI;

   R=I-PAI(+1);

   D_GDP_NW=Y-Y(-1)+ZGDP;

   DPQ_P_NW=paiss+PAI;

   RN3M_NW=iss+I;

   ZI=rhoi*ZI(-1)+sigi*EI;
   
   ZPAI=rhopai*ZPAI(-1)+sigpai*EPAI;
   
   ZY=rhoy*ZY(-1)+sigy*EY;
   
   ZGDP=(1-rhogdp)*gyss+rhogdp*ZGDP(-1)+siggdp*EGDP;


parameterization
% not estimated
	gyss   		 ,0 	      ;						 
	iss    		 ,0 	      ;						 
	paiss  		 ,0 	      ;
	beta_lag 	 ,0.5000      ;						 
	beta_lead	 ,0.4000	  ;						 
	beta_r  	 ,0.9000	  ;						 
	gam_lag 	 ,0.6000	  ;						 
	gam_y   	 ,0.5000	  ;						 
	lamb_lag	 ,0.8000	  ;						 
	lamb_lead    ,0.1000	  ;						 
	lamb_y  	 ,0.3000	  ;						 
	rhogdp 		 ,0.5000	  ;						 
	rhoi   		 ,0.5000	  ;						 
	rhopai 		 ,0.5000	  ;						 
	rhoy   		 ,0.5000      ;						 
	siggdp 		 ,0.5000	  ;% 0.3138  ,14.1339						 
	sigi   		 ,0.5000	  ;						 
	sigpai 		 ,0.5000	  ;						 
	sigy   		 ,0.5000	  ;

	