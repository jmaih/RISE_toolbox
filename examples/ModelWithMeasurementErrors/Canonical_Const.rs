endogenous	DPQ_P_NW "Inflation", Y, ZGDP, ZI, ZPAI, ZY, D_GDP_NW "Growth", I, PAI, R, RN3M_NW "Fed Funds Rate"

exogenous EGDP "output shock",EI "monetary policy shock",EPAI "Cost push shock",EY "IS shock"

parameters beta_lag "\beta_{lag}", beta_lead "\beta_{lead}", beta_r "\beta_{r}",
gam_lag "\gamma_{lag}", gam_y "\gamma_{y}", gyss, iss, lamb_lag "\lambda_{lag}",
lamb_lead "\lambda_{lead}", lamb_y "\lambda_{y}", paiss,
rhogdp "\rho_{gdp}", rhoi "\rho_{i}", rhopai "\rho_{\pi}", rhoy "\rho_{y}",
siggdp "\sigma_{gdp}", sigi "\sigma_{i}", sigpai "\sigma_{\pi}", sigy "\sigma_{y}"

% measurement errors are declared as parameters: The corresponding variables
% have to be declared as observables!!!
%---------------------------------------------------------------------------
parameters stderr_DPQ_P_NW, stderr_D_GDP_NW, stderr_RN3M_NW 


model(linear)
	
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

   
observables DPQ_P_NW, D_GDP_NW, RN3M_NW;


parameterization
	gyss   		 ,0 	      ;						 
	iss    		 ,0 	      ;						 
	paiss  		 ,0 	      ;
	beta_lag 	 ,0.5000      ,0.1     ,.9000,gam_pdf,.90;						 
	beta_lead	 ,0.4000	  ,0.1     ,.9000,gam_pdf,.90;						 
	beta_r  	 ,0.9000	  ,0.1     ,3.000,gam_pdf,.90;						 
	gam_lag 	 ,0.6000	  ,0.1     ,.8000,beta_pdf,.90;						 
	gam_y   	 ,0.5000	  ,0.1     ,.9000,gam_pdf,.90;						 
	lamb_lag	 ,0.8000	  ,0.1     ,.9000,gam_pdf,.90;						 
	lamb_lead    ,0.1000	  ,0.1     ,.9000,gam_pdf,.90;						 
	lamb_y  	 ,0.3000	  ,0.1     ,2.000,gam_pdf,.90;						 
	rhogdp 		 ,0.5000	  ,.1      ,.8000,beta_pdf,.90;						 
	rhoi   		 ,0.5000	  ,.1      ,.8000,beta_pdf,.90;						 
	rhopai 		 ,0.5000	  ,.1      ,.8000,beta_pdf,.90;						 
	rhoy   		 ,0.5000      ,.1      ,.8000,beta_pdf,.90;						 
	siggdp 		 ,0.5000	  ,0.0050  ,1.000,inv_gamma_pdf,.90;						 
	sigi   		 ,0.5000	  ,0.0050  ,1.000,inv_gamma_pdf,.90;						 
	sigpai 		 ,0.5000	  ,0.0050  ,1.000,inv_gamma_pdf,.90;						 
	sigy   		 ,0.5000	  ,0.0050  ,1.000,inv_gamma_pdf,.90;
%  Measurement errors
	stderr_DPQ_P_NW   		 ,0.5000	  ,0.0050  ,1.000,inv_gamma_pdf,.90;
	stderr_D_GDP_NW   		 ,0.5000	  ,0.0050  ,1.000,inv_gamma_pdf,.90;
	stderr_RN3M_NW   		 ,0.5000	  ,0.0050  ,1.000,inv_gamma_pdf,.90;

	