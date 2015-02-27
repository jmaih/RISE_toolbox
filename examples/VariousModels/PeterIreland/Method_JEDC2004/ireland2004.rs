% Peter N. Ireland: "A Method for Taking Models to the Data" 
% Journal of Economic Dynamics & Control 28 (2004) 1205 – 1226

endogenous A Y K I C H LY LC LH
U_Y U_C U_H XI_Y XI_C XI_H

exogenous EPS EPS_Y EPS_C EPS_H	TREND

parameters theta a eta delta gam beta rho sigma
dyy dyc dyh                                            
dcy dcc dch
dhy dhc dhh
vyy
vcy vcc
vhy vhc vhh

observables LY LC LH TREND

model

	Y=A*K{-1}^theta*H^(1-theta);

	log(A)=(1-rho)*log(a)+rho*log(A{-1})+sigma*EPS;

	Y=C+I;

	eta*K=(1-delta)*K{-1}+I;

	gam*C*H=(1-theta)*Y;

	eta/C=beta*1/C{+1}*(theta*Y{+1}/K+1-delta);

	% Measurement equations					                                 
	%----------------------					                                 
	LY = log(Y) + log(eta)*TREND + U_Y; 
	%  ythat = log(yt) - log(eta)*trend' - log(yss); 

	LC = log(C) + log(eta)*TREND + U_C;
	%  cthat = log(ct) - log(eta)*trend' - log(css); 

	LH=log(H)+ U_H;
	%  hthat = log(ht) - log(hss);                   

	% VAR model for measurement errors
	%-----------------------------------
	U_Y = dyy*U_Y{-1}+ dyc*U_C{-1}+ dyh*U_H{-1}+XI_Y; 
   
	U_C = dcy*U_Y{-1}+ dcc*U_C{-1}+ dch*U_H{-1}+XI_C; 
   
	U_H = dhy*U_Y{-1}+ dhc*U_C{-1}+ dhh*U_H{-1}+XI_H; 

	% correlated shocks in Cholesky form
	%------------------------------------
	XI_Y = vyy*EPS_Y;
    
	XI_C = vcy*EPS_Y+ vcc*EPS_C;
    
	XI_H = vhy*EPS_Y+ vhc*EPS_C+ vhh*EPS_H;
    
steady_state_model
	xx_ssmdef_1=(eta/beta-1+delta);
	xx_ssmdef_2=(eta-1+delta);
	A=a;
	Y=A^(1/(1-theta))*(theta/xx_ssmdef_1)^(theta/(1-theta))*(1-theta)/gam*1/(1-theta*xx_ssmdef_2/xx_ssmdef_1);
	K=theta/xx_ssmdef_1*Y;
	I=theta*xx_ssmdef_2/xx_ssmdef_1*Y;
	C=(1-theta*xx_ssmdef_2/xx_ssmdef_1)*Y;
	H=(1-theta)/gam*1/(1-theta*xx_ssmdef_2/xx_ssmdef_1);
	LY = log(Y); 
	LC = log(C); 
	LH = log(H);

parameterization
	delta, 0.025;
	beta,  0.99;
	theta, 0.2292, 0, 1;
	a, 5.1847, 0, 10;
	eta, 1.0051, 1, 3;
	gam, 0.0045, 0, 3;
	rho, 0.9987, 0, 1;
	sigma, 0.0056, 0, 3;
	dyy ,1.3655-1, -2, 2 ;
	dyc ,0.3898*0, -2, 2  ;
	dyh ,0.4930*0, -2, 2  ;
	dcy ,0.1380*0, -2, 2  ;
	dcc ,0.9690, -2, 2  ;
	dch ,0.1046*0, -2, 2  ;
	dhy ,0.7153*0, -2, 2  ;
	dhc ,0.4605*0, -2, 2  ;
	dhh ,0.2219, -2, 2  ;
	vyy , 0.007, 0, 3;
	vcy , 0.0042701, -2, 2  ;
	vcc , 0.0054202, 0, 3;
	vhy , 0.0012898, -2, 2  ;
	vhc , 0.0012645, -2, 2  ;
	vhh , 0.00012545, 0, 3;
	
%V=[0.0070^2 0.00002989 0.00000903
%0.00002989 0.0069^2 0.00001237
%0.00000903 0.00001237 0.0018^2]
% not positive definnite: nearest
% covariance had to be used
