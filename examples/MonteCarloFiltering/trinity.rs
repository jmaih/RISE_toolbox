endogenous	PIE, GROWTH, I, PAI, RS, Y

exogenous EI, EPAI, EY 

parameters pc_beta, d_y_lead, d_r,
mp_i_lag, mp_pai, mp_y, gyss, iss,
pc_pai_lead, pc_y, paiss, sigi, sigpai, sigy

model
	% Phillips curve
   PAI=pc_beta*pc_pai_lead*PAI(+1)+(1-pc_pai_lead)*PAI(-1)+pc_y*Y(-1)+sigpai*EPAI;

	% Demand curve
   Y=(1-d_y_lead)*Y(-1)+d_y_lead*Y(+1)-d_r*(I-PAI(+1))+sigy*EY;

	% Monetary policy reaction
   I=mp_i_lag*I(-1)+(1-mp_i_lag)*(mp_pai*PAI+mp_y*Y)+sigi*EI;

	% Measurement equations
   GROWTH=gyss+400*(Y-Y(-1));	% +ZGDP
   
   PIE=paiss+400*PAI;
   
   RS=iss+400*I;

observables PIE, GROWTH, RS


steady_state_model
%%	ZGDP=gyss;
	GROWTH=gyss;
	PIE=paiss;
	RS=iss;

parameterization
% Not estimated
	pc_beta       ,0.9900;				 
% Controlled by the constant chain
	gyss          ,1.0000     ,0       ,5  ;				 
	iss 	      ,5.0000     ,0       ,30 ;						 
	paiss         ,5.0000     ,0       ,30 ;
	pc_pai_lead   ,0.0100	  ,0	   ,.99;						 
	pc_y   	      ,0.3000	  ,0	   ,1.5;
	d_y_lead      ,0.7000	  ,0	   ,.99;						 
	d_r   	      ,0.0600	  ,0	   ,3.0;						 
	mp_i_lag  	  ,0.7000	  ,0	   ,.99;
	mp_pai    	  ,2.5000	  ,1	   ,3.0;
	mp_y    	  ,0.5000	  ,0	   ,.99;						 
	sigi 	      ,0.0500	  ,0.0001  ,1.0;						 
	sigpai 	      ,0.0500	  ,0.0001  ,1.0;						 
	sigy 	      ,0.0500	  ,0.0001  ,1.0;
