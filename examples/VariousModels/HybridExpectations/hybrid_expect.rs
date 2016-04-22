endogenous	Y, ZI, ZPAI, ZY, I, PAI, R

exogenous EI "monetary policy shock",EPAI "Cost push shock",EY "IS shock"

parameters beta_lag "$\beta_{lag}$", beta_lead "$\beta_{lead}$", beta_r "$\beta_{r}$",
gam_lag "$\gamma_{lag}$", gam_y "$\gamma_{y}$", lamb_lag "$\lambda_{lag}$",
lamb_lead "$\lambda_{lead}$", lamb_y "$\lambda_{y}$",
rhoi "$\rho_{i}$", rhopai "$\rho_{\pi}$", rhoy "$\rho_{y}$",
siggdp "$\sigma_{gdp}$", sigi "$\sigma_{i}$", sigpai "$\sigma_{\pi}$", sigy "$\sigma_{y}$"

parameters hbe_lambda hbe_w

model

   Y=beta_lag*Y(-1)+beta_lead*Y(+1)-beta_r*R(-1)+ZY;

   PAI=lamb_lag*PAI(-1)+lamb_lead*PAI(+1)+lamb_y*Y(-1)+ZPAI;

   I=gam_lag*I(-1)+(1-gam_lag)*(PAI(+4)+gam_y*Y)+ZI;

   R=I-PAI(+1);

   ZI=rhoi*ZI(-1)+sigi*EI;
   
   ZPAI=rhopai*ZPAI(-1)+sigpai*EPAI;
   
   ZY=rhoy*ZY(-1)+sigy*EY;

   
parameterization;
	beta_lag 	 ,0.0000;
	beta_lead	 ,0.9900;
	beta_r  	 ,0.1000;
	gam_lag 	 ,0.6000;
	gam_y   	 ,0.5000;
	lamb_lag	 ,0.0000;
	lamb_lead    ,0.9900;
	lamb_y  	 ,0.0500;
	rhoi   		 ,0.0000;
	rhopai 		 ,0.0000;
	rhoy   		 ,0.0000;
	siggdp 		 ,0.5000;					 
	sigi   		 ,0.5000; 
	sigpai 		 ,0.5000; 
	sigy   		 ,0.5000;
	hbe_lambda   ,0.3000; % weight on current relative to past in backward-looking expectations
	hbe_w        ,0.3000; % weight on backward-looking expectations in "total" expectations

	