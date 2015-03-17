endogenous C I AG H AL G_H A ZG_C LAMBDA_C H_C K_C I_C Z_C XI_C ZL_C G_C ZG_I LAMBDA_I H_I K_I I_I Z_I XI_I ZL_I G_I
LG_C LG_I LG_H

observables LG_C LG_I LG_H

exogenous ELA EAG EL_C EG_C EL_I EG_I

parameters gam beta	rhola sigla	agtr rhoga sigga  phik_C_tr phih_C_tr theta_C eta_C	kappa_C	delta_C	rhol_C
	sigl_C zg_C_tr rhog_C sigg_C theta_I eta_I kappa_I delta_I rhol_I sigl_I zg_I_tr rhog_I sigg_I
	
model
	# phik_C = 100*abs(phik_C_tr);
	# phik_I = phik_C;
	# phih_C = 100*abs(phih_C_tr);
	# phih_I = phih_C;
	# ag = 1 + agtr;
	# zg_C = 1 + abs(zg_C_tr);
	# zg_I = 1 + abs(zg_I_tr);
	
	1/(C-gam*C{-1}/(AG{-1}*ZG_I{-1}^theta_C*ZG_C{-1}^(1-theta_C)))-
		beta*gam*1/((AG*ZG_I^theta_C*ZG_C^(1-theta_C))*C{+1}-gam*C)=LAMBDA_C;
		
	1/A=(1-theta_C)*LAMBDA_C*C/H_C-phih_C*LAMBDA_C*(AG{-1}*H_C/H_C{-1}-eta_C)*
		AG{-1}/H_C{-1}*(1-phik_C/2*(I_C/K_C{-1}-kappa_C)^2)*K_C{-1}^theta_C*(Z_C*H_C)^(1-theta_C)+
		beta*phih_C*LAMBDA_C{+1}*(AG*H_C{+1}/H_C-eta_C)*AG*H_C{+1}/H_C*1/H_C*(1-phik_C/2*
		(I_C{+1}/K_C-kappa_C)^2)*K_C^theta_C*(Z_C{+1}*H_C{+1})^(1-theta_C);
		
	1/A=(1-theta_I)*LAMBDA_I*I/H_I-phih_I*LAMBDA_I*(AG{-1}*H_I/H_I{-1}-eta_I)*
		AG{-1}/H_I{-1}*(1-phik_I/2*(I_I/K_I{-1}-kappa_I)^2)*K_I{-1}^theta_I*(Z_I*H_I)^(1-theta_I)+
		beta*phih_I*LAMBDA_I{+1}*(AG*H_I{+1}/H_I-eta_I)*AG*H_I{+1}/H_I*1/H_I*(1-phik_I/2*
		(I_I{+1}/K_I-kappa_I)^2)*K_I^theta_I*(Z_I{+1}*H_I{+1})^(1-theta_I);
		
	XI_C=LAMBDA_I+phik_C*LAMBDA_C*(1-phih_C/2*(AG{-1}*H_C/H_C{-1}-eta_C)^2)*
			(I_C/K_C{-1}-kappa_C)*K_C{-1}^(theta_C-1)*(Z_C*H_C)^(1-theta_C);
			
	XI_I=LAMBDA_I*(1+phik_I*(1-phih_I/2*(AG{-1}*H_I/H_I{-1}-eta_I)^2)*
			(I_I/K_I{-1}-kappa_I)*K_I{-1}^(theta_I-1)*(Z_I*H_I)^(1-theta_I));
		
	AG*Z_I*XI_C=beta*(1-delta_C)*XI_C{+1}+beta*theta_C*LAMBDA_C{+1}*C{+1}/K_C+beta*phik_C*LAMBDA_C{+1}*
		(1-phih_C/2*(AG*H_C{+1}/H_C-eta_C)^2)*
			(I_C{+1}/K_C-kappa_C)*I_C{+1}/K_C*K_C^(theta_C-1)*(Z_C{+1}*H_C{+1})^(1-theta_C);
			
	AG*Z_I*XI_I=beta*(1-delta_I)*XI_I{+1}+beta*theta_I*LAMBDA_I{+1}*I{+1}/K_I+beta*phik_I*LAMBDA_I{+1}*
		(1-phih_I/2*(AG*H_I{+1}/H_I-eta_I)^2)*
			(I_I{+1}/K_I-kappa_I)*I_I{+1}/K_I*K_I^(theta_I-1)*(Z_I{+1}*H_I{+1})^(1-theta_I);
			
	C = (1-phih_C/2*(AG{-1}*H_C/H_C{-1}-eta_C)^2)*(1-phik_C/2*(I_C/K_C{-1}-kappa_C)^2)*K_C{-1}^theta_C*(Z_C*H_C)^(1-theta_C);
	
	I = (1-phih_I/2*(AG{-1}*H_I/H_I{-1}-eta_I)^2)*(1-phik_I/2*(I_I/K_I{-1}-kappa_I)^2)*K_I{-1}^theta_I*(Z_I*H_I)^(1-theta_I);
	
	(1-delta_C)*K_C{-1}+I_C=AG*ZG_I*K_C;
	
	(1-delta_I)*K_I{-1}+I_I=AG*ZG_I*K_I;
	
	I=I_C+I_I;
	
	H=H_C+H_I;
	
	log(A)=log(AL)+log(AG);
	
	log(AL)=rhola*log(AL{-1})+sigla*ELA;
	
	log(AG/ag)=rhoga*log(AG{-1}/ag)+sigga*EAG;
	
	log(Z_C)=log(ZL_C)+log(ZG_C);
	
	log(ZL_C)=rhol_C*log(ZL_C{-1})+sigl_C*EL_C;
	
	log(ZG_C/zg_C)=rhog_C*log(ZG_C{-1}/zg_C)+sigg_C*EG_C;
	
	log(Z_I)=log(ZL_I)+log(ZG_I);
	
	log(ZL_I)=rhol_I*log(ZL_I{-1})+sigl_I*EL_I;
	
	log(ZG_I/zg_I)=rhog_I*log(ZG_I{-1}/zg_I)+sigg_I*EG_I;
	
	G_C=AG{-1}*ZG_I{-1}^theta_C*ZG_C{-1}^(1-theta_C)*C/C{-1};
	
	G_I=AG{-1}*ZG_I{-1}*I/I{-1};
	
	G_H=AG{-1}*H/H{-1};
	
	LG_C=log(G_C);
	
	LG_I=log(G_I);
	
	LG_H=log(G_H);
