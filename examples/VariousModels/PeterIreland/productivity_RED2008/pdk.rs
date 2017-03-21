endogenous C "Consumption" I "Total investment" AG "Growth-rate-stationary pref. process" H "Total labor supply"
AL "level-stationary pref. process" G_H "Growth rate of hours" A "Preference shock"

@#for ii in Sektors		
	ZG_@{ii} "Growth-rate-stationary techn. process @{ii}"
	LAMBDA_@{ii} "Multiplier on production @{ii}"
	H_@{ii} "Labor supply in sector @{ii}"
	K_@{ii} "Capital stock in sector @{ii}"
	I_@{ii}	"Investment goods in sector @{ii}"
	Z_@{ii} "Technology shock in sector @{ii}"
	XI_@{ii} "Multiplier on capital accumulation @{ii}"
	ZL_@{ii} "Level-stationary techn. process @{ii}"
	G_@{ii} "Growth rate of @{ii}"
@#end

LG_C LG_I LG_H

observables LG_C LG_I LG_H

exogenous ELA EAG
@#for ii in Sektors		
	EL_@{ii} EG_@{ii}
@#end

parameters gam beta	rhola sigla	agtr rhoga sigga  phik_C_tr phih_C_tr

@#for ii in Sektors
	theta_@{ii} eta_@{ii}	kappa_@{ii}	delta_@{ii}	rhol_@{ii}
	sigl_@{ii} zg_@{ii}_tr	rhog_@{ii}	sigg_@{ii}
@#end

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

@#for ii in Sektors		
	1/A=(1-theta_@{ii})*LAMBDA_@{ii}*@{ii}/H_@{ii}-phih_@{ii}*LAMBDA_@{ii}*(AG{-1}*H_@{ii}/H_@{ii}{-1}-eta_@{ii})*
		AG{-1}/H_@{ii}{-1}*(1-phik_@{ii}/2*(I_@{ii}/K_@{ii}{-1}-kappa_@{ii})^2)*K_@{ii}{-1}^theta_@{ii}*(Z_@{ii}*H_@{ii})^(1-theta_@{ii})+
		beta*phih_@{ii}*LAMBDA_@{ii}{+1}*(AG*H_@{ii}{+1}/H_@{ii}-eta_@{ii})*AG*H_@{ii}{+1}/H_@{ii}*1/H_@{ii}*(1-phik_@{ii}/2*
		(I_@{ii}{+1}/K_@{ii}-kappa_@{ii})^2)*K_@{ii}^theta_@{ii}*(Z_@{ii}{+1}*H_@{ii}{+1})^(1-theta_@{ii});

	XI_@{ii}=LAMBDA_I
		@#if strcmp('@{ii}','I')
			*(1
		@#end
		+phik_@{ii}*
		@#if strcmp('@{ii}','C')
			LAMBDA_@{ii}*
		@#end
		(1-phih_@{ii}/2*(AG{-1}*H_@{ii}/H_@{ii}{-1}-eta_@{ii})^2)*
			(I_@{ii}/K_@{ii}{-1}-kappa_@{ii})*K_@{ii}{-1}^(theta_@{ii}-1)*(Z_@{ii}*H_@{ii})^(1-theta_@{ii})
		@#if strcmp('@{ii}','I')
			)
		@#end
		;

	AG*Z_I*XI_@{ii}=beta*(1-delta_@{ii})*XI_@{ii}{+1}+beta*theta_@{ii}*LAMBDA_@{ii}{+1}*@{ii}{+1}/K_@{ii}+beta*phik_@{ii}*LAMBDA_@{ii}{+1}*
		(1-phih_@{ii}/2*(AG*H_@{ii}{+1}/H_@{ii}-eta_@{ii})^2)*
			(I_@{ii}{+1}/K_@{ii}-kappa_@{ii})*I_@{ii}{+1}/K_@{ii}*K_@{ii}^(theta_@{ii}-1)*(Z_@{ii}{+1}*H_@{ii}{+1})^(1-theta_@{ii});

	@{ii} = (1-phih_@{ii}/2*(AG{-1}*H_@{ii}/H_@{ii}{-1}-eta_@{ii})^2)*(1-phik_@{ii}/2*(I_@{ii}/K_@{ii}{-1}-kappa_@{ii})^2)*K_@{ii}{-1}^theta_@{ii}*(Z_@{ii}*H_@{ii})^(1-theta_@{ii});

	(1-delta_@{ii})*K_@{ii}{-1}+I_@{ii}=AG*ZG_I*K_@{ii};
	
	log(Z_@{ii})=log(ZL_@{ii})+log(ZG_@{ii});
	
	log(ZL_@{ii})=rhol_@{ii}*log(ZL_@{ii}{-1})+sigl_@{ii}*EL_@{ii};

	log(ZG_@{ii}/zg_@{ii})=rhog_@{ii}*log(ZG_@{ii}{-1}/zg_@{ii})+sigg_@{ii}*EG_@{ii};
@#end

	I=I_C+I_I;

	H=H_C+H_I;

	log(A)=log(AL)+log(AG);

	log(AL)=rhola*log(AL{-1})+sigla*ELA;

	log(AG/ag)=rhoga*log(AG{-1}/ag)+sigga*EAG;

	G_C=AG{-1}*ZG_I{-1}^theta_C*ZG_C{-1}^(1-theta_C)*C/C{-1};

	G_I=AG{-1}*ZG_I{-1}*I/I{-1};

	G_H=AG{-1}*H/H{-1};

	% measurement equations
	%-----------------------

	LG_C=log(G_C);
	
	LG_I=log(G_I);
	
	LG_H=log(G_H);

%	LG_C=log(steady_state(G_C))*(1+log(G_C/steady_state(G_C)));
%	
%	LG_I=log(steady_state(G_I))*(1+log(G_I/steady_state(G_I)));
%	
%	LG_H=log(steady_state(G_H))*(1+log(G_H/steady_state(G_H)));



