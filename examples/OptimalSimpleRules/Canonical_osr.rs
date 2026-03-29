@endogenous	Y, I, PAI, R, DI

@exogenous EI, EPAI, EY

@parameters beta_lag, beta_lead, beta_r, gam_lag, gam_y, lamb_lag,
lamb_lead, lamb_y, sigi, sigpai, sigy

@model

   Y=beta_lag*Y(-1)+beta_lead*Y(+1)-beta_r*R(-1)+sigy*EY;

   PAI=lamb_lag*PAI(-1)+lamb_lead*PAI(+1)+lamb_y*Y(-1)+sigpai*EPAI;

   % Taylor rule
   I=gam_lag*I(-1)+(1-gam_lag)*(PAI(+4)+gam_y*Y)+sigi*EI;
				  
   R=I-PAI(+1);
   % interest rate smoothing
   DI=I-I(-1);
