% Peter N. Ireland(2013): "Stochastic Growth In the United States and Euro Area"
% Journal of the European Economic Association vol. 11 number 1, pages1:24
%
% Model 1: Nonstationary and Cointegrated
%	- Neutral Technology Shocks
%	- Investment-Specific Technology Shocks
% Model 2 = Model 1 + Nonstationary and Cointegrated
%	- Preference Shocks
% Two-Country Model with Capital and Labor Adjustment Costs
%
% start date: 1970Q1

endogenous
@#for cc in CountryNames
	C_@{cc}, "Consumption(@{cc})" I_@{cc}, "Investment(@{cc})" G_@{cc},  "Gov. Spending process(@{cc})" N_@{cc}, "Net Exports(@{cc})"
	CTILDE_@{cc}, "Total consumption(@{cc})" ITILDE_@{cc}, "Total investment(@{cc})" K_@{cc}, "Capital(@{cc})"
	D_@{cc}, "Bonds(@{cc})" LAMBDA_@{cc} "Lagrange mult. on budget(@{cc})" XI_@{cc}, "Lagrange mult. on capital accum.(@{cc})"
	TAU_@{cc}, "Lump-sum taxes(@{cc})" A_@{cc}, "Demand of home interm. goods(@{cc})" 
	B_@{cc}, "Demand of Foreign interm. goods(@{cc})" P_@{cc}, "Price of cons. goods(@{cc})" 
	X_@{cc}, "Price of Inv. goods(@{cc})" W_@{cc}, "Nom. wages(@{cc})" 
	Q_@{cc}, "Nom. rental rate of capital(@{cc})" M_@{cc}, "Preference shock process(@{cc})" 
	Z_@{cc}, "Neutral tech. shock process(@{cc})", V_@{cc} "Inv.-specific tech. shock process(@{cc})"
	L_@{cc}, "Hours worked(@{cc})" R_GC_@{cc} "Gov. spending-Consumption ratio(@{cc})" 
@#end
Y_A, "Home intermediate goods" Y_B "Foreign intermediate goods", R "Interest rate",
P_A, "real price of home interm. goods" P_B, "real price of foreign interm. goods" RER, "real exch. rate (pF)"
TOT, "Terms of trade(pB/pA)" Z_HF, "Neutral tech. ratio(H/F)" V_HF, "Inv.-specific tech. ratio(H/F)" G_CH, "Consumption growth"
G_IH, "Inv. growth" R_CFH, "Foreign-Home consumption ratio" R_IFH "Foreign/Home Invest. ratio",
@#if ncp_shocks
	M_HF "Preferences shocks ratio (H/F)", G_LH "Growth of Hours worked ratio (H)", R_LFH "Hours worked ratio (F/H)"
	G_LH_HAT R_LFH_HAT
@#end
G_CH_HAT, G_IH_HAT, R_GC_H_HAT, L_H_HAT, R_CFH_HAT, R_IFH_HAT, R_GC_F_HAT, L_F_HAT

observables G_CH_HAT, G_IH_HAT, R_GC_H_HAT, R_CFH_HAT, R_IFH_HAT, R_GC_F_HAT,
@#if ncp_shocks
	G_LH_HAT R_LFH_HAT
@#else
	 L_H_HAT, L_F_HAT
@#end

exogenous
@#for cc in CountryNames
	EM_@{cc} "Preference(@{cc})" EV_@{cc} "Inv.-specific tech(@{cc})" EZ_@{cc} "Neutral tech(@{cc})" EG_@{cc} "Gov. spending(@{cc})"
@#end

parameters
alpha, phi_d philtr phiktr mu gam delta theta omega z_hf_ss v_hf_ss beta kappazhtr
kappavhtr
@#if ncp_shocks
m_hf_ss
@#end

@#for cc in CountryNames
	eta_@{cc},
	rho_m_@{cc}@{cc}, rho_z_@{cc}@{cc}, rho_v_@{cc}@{cc}, rho_g_@{cc}@{cc}, 
	@#if strcmp('@{cc}','H')
		 rho_m_@{cc}F,rho_z_@{cc}F, rho_v_@{cc}F, rho_g_@{cc}F, 
	@#else
		rho_m_@{cc}H, rho_z_@{cc}H, rho_v_@{cc}H, rho_g_@{cc}H
		kappa_z_@{cc}, kappa_v_@{cc}
	@#end
	mss_@{cc}, zss_@{cc}, vss_@{cc}, gss_@{cc}
	sig_g_@{cc}	sig_v_@{cc}	sig_z_@{cc} sig_m_@{cc}
	@#if ncp_shocks	&& strcmp('@{cc}','F')
		kappa_m_@{cc} kappamhtr
	@#end
@#end

model
	# phi_l = 100*abs(philtr);
	# phi_k = 100*abs(phiktr);
	# kappa_z_H = -abs(kappazhtr);
	# kappa_v_H = -abs(kappavhtr);
	@#if ncp_shocks
		# kappa_m_H = -abs(kappamhtr);
	@#end
@#for cc in CountryNames
	W_@{cc}*L_@{cc}+Q_@{cc}*K_@{cc}{-1}+D_@{cc}{-1}
	@#if strcmp('@{cc}','F')
		/P_F
	@#end
	= C_@{cc}+X_@{cc}*I_@{cc}+TAU_@{cc}+
	@#if ncp_shocks
		M_@{cc}*
	@#end
	V_@{cc}^(alpha/(1-alpha))*Z_@{cc}*D_@{cc}/(
	@#if strcmp('@{cc}','F')
		P_F*
	@#end
	R)+phi_d/2*
	@#if ncp_shocks
		M_@{cc}*
	@#end
	V_@{cc}^(alpha/(1-alpha))*Z_@{cc}*D_@{cc}^2+phi_l/2*(
	@#if ncp_shocks
		M_@{cc}{-1}*
	@#end
	L_@{cc}/L_@{cc}{-1}-
	@#if ncp_shocks
		steady_state(M_@{cc})
	@#else
		1
	@#end
	)^2*L_@{cc}{-1}
	@#if ncp_shocks
		/M_@{cc}{-1}
	@#end % @#else creates wrong error message...
	;

	(1-delta)*K_@{cc}{-1}+I_@{cc}-phi_k/2*(I_@{cc}/K_@{cc}{-1}-eta_@{cc})^2*K_@{cc}{-1}=
	@#if ncp_shocks
		M_@{cc}*
	@#end
	V_@{cc}^(1/(1-alpha))*Z_@{cc}*K_@{cc};

	mu*(C_@{cc}^mu*(1-L_@{cc}/M_@{cc})^(1-mu))^(1-gam)=LAMBDA_@{cc}*C_@{cc};

	(1-mu)*(C_@{cc}^mu*(1-L_@{cc}/M_@{cc})^(1-mu))^(1-gam)/(M_@{cc}*(1-L_@{cc}/M_@{cc}))=
		LAMBDA_@{cc}*(W_@{cc}-phi_l*(
		@#if ncp_shocks
			M_@{cc}{-1}*
		@#end
		L_@{cc}/L_@{cc}{-1}-
		@#if ncp_shocks
			steady_state(M_@{cc})
		@#else
			1
		@#end
		))+
		beta*phi_l*LAMBDA_@{cc}{+1}
		@#if ncp_shocks
			/M_@{cc}
		@#end
		*(
		@#if ncp_shocks
			M_@{cc}*
		@#end
		V_@{cc}^(alpha/(1-alpha))*Z_@{cc})^(mu*(1-gam))*((
		@#if ncp_shocks
			M_@{cc}*
		@#end
		L_@{cc}{+1}/L_@{cc}-
		@#if ncp_shocks
			steady_state(M_@{cc})
		@#else
			1
		@#end
		)*
		@#if ncp_shocks
			M_@{cc}*
		@#end
		L_@{cc}{+1}/L_@{cc}
		-1/2*(
		@#if ncp_shocks
			M_@{cc}*
		@#end
		L_@{cc}{+1}/L_@{cc}-
		@#if ncp_shocks
			steady_state(M_@{cc})
		@#else
			1
		@#end
		)^2);

	LAMBDA_@{cc}*X_@{cc} = XI_@{cc}*(1-phi_k*(I_@{cc}/K_@{cc}{-1}-eta_@{cc}));

	(
	@#if ncp_shocks
		M_@{cc}*
	@#end
	V_@{cc}^(alpha/(1-alpha))*Z_@{cc})^(1-mu*(1-gam))*V_@{cc}*XI_@{cc} = beta*LAMBDA_@{cc}{+1}*Q_@{cc}{+1} +
		beta*XI_@{cc}{+1}*(1-delta+phi_k*(I_@{cc}{+1}/K_@{cc}-eta_@{cc})*I_@{cc}{+1}/K_@{cc}-phi_k/2*(I_@{cc}{+1}/K_@{cc}-eta_@{cc})^2);

	(
	@#if ncp_shocks
		M_@{cc}*
	@#end
	V_@{cc}^(alpha/(1-alpha))*Z_@{cc})^(1-mu*(1-gam))*LAMBDA_@{cc}*(1/(
	@#if strcmp('@{cc}','F')
		P_F*
	@#end
	R)+phi_d*D_@{cc})=beta*LAMBDA_@{cc}{+1}
	@#if strcmp('@{cc}','F')
		/P_F{+1}
	@#end
	;

	K_@{cc}{-1}^alpha*(Z_@{cc}*L_@{cc})^(1-alpha)=
	@#if strcmp('@{cc}','H')
		Y_A
	@#else
		Y_B
	@#end
	;
		
	alpha*	
	@#if strcmp('@{cc}','H')
		P_A*Y_A
	@#else
		P_B*Y_B
	@#end
		=
	@#if strcmp('@{cc}','F')
		P_F*
	@#end
		Q_@{cc}*K_@{cc}{-1};
		
	(1-alpha)*	
	@#if strcmp('@{cc}','H')
		P_A*Y_A
	@#else
		P_B*Y_B
	@#end
		=
	@#if strcmp('@{cc}','F')
		P_F*
	@#end
		W_@{cc}*L_@{cc};

	(
	@#if strcmp('@{cc}','H')
		(1-omega)
	@#else
		omega
	@#end
	^(1/theta)*A_@{cc}^((theta-1)/theta)+
	@#if strcmp('@{cc}','H')
		omega
	@#else
		(1-omega)
	@#end
	^(1/theta)*B_@{cc}^((theta-1)/theta))^(theta/(theta-1))=CTILDE_@{cc}+1/V_@{cc}*ITILDE_@{cc};

	X_@{cc} = 1/V_@{cc};

		P_A
	@#if strcmp('@{cc}','F')
		/P_F
	@#end
		= (
	@#if strcmp('@{cc}','H')
		(1-omega)
	@#else
		omega
	@#end
		^(1/theta)*A_@{cc}^((theta-1)/theta)+
	@#if strcmp('@{cc}','H')
		omega
	@#else
		(1-omega)
	@#end
		^(1/theta)*B_@{cc}^((theta-1)/theta))^(1/(theta-1))*
	@#if  strcmp('@{cc}','H')
		(1-omega)
	@#else
		omega
	@#end
		^(1/theta)*A_@{cc}^(-1/theta);

		P_B
	@#if  strcmp('@{cc}','F')
		/P_F
	@#end
		= (
	@#if  strcmp('@{cc}','H')
		(1-omega)
	@#else
		omega
	@#end
		^(1/theta)*A_@{cc}^((theta-1)/theta)+
	@#if  strcmp('@{cc}','H')
		omega
	@#else
		(1-omega)
	@#end
		^(1/theta)*B_@{cc}^((theta-1)/theta))^(1/(theta-1))*
	@#if  strcmp('@{cc}','H')
		omega
	@#else
		(1-omega)
	@#end
		^(1/theta)*B_@{cc}^(-1/theta);

	TAU_@{cc} = G_@{cc};

	CTILDE_@{cc} = C_@{cc}+phi_d/2*
		@#if ncp_shocks
			M_@{cc}*
		@#end
		V_@{cc}^(alpha/(1-alpha))*Z_@{cc}*D_@{cc}^2+
		phi_l/2*(
		@#if ncp_shocks
			M_@{cc}{-1}*
		@#end
		L_@{cc}/L_@{cc}{-1}-
		@#if ncp_shocks
			steady_state(M_@{cc})
		@#else
			1
		@#end
		)^2*L_@{cc}{-1}
		@#if ncp_shocks
			/M_@{cc}{-1}
		@#end
		+G_@{cc};

	ITILDE_@{cc} = I_@{cc};

	log(M_@{cc}/mss_@{cc})=rho_m_@{cc}@{cc}*log(M_@{cc}{-1}/mss_@{cc})+ 	
	@#if  strcmp('@{cc}','H')
		rho_m_HF*log(M_F{-1}/mss_F)
	@#else
		rho_m_FH*log(M_H{-1}/mss_H)
	@#end
	@#if ncp_shocks
		+kappa_m_@{cc}*log(M_HF{-1}/m_hf_ss)
	@#end
	+sig_m_@{cc}*EM_@{cc};

	log(Z_@{cc}/zss_@{cc})=rho_z_@{cc}@{cc}*log(Z_@{cc}{-1}/zss_@{cc})	
	@#if  strcmp('@{cc}','H')
		+ rho_z_HF*log(Z_F{-1}/zss_F)
	@#else
		+ rho_z_FH*log(Z_H{-1}/zss_H)
	@#end
	+kappa_z_@{cc}*log(Z_HF/z_hf_ss)+sig_z_@{cc}*EZ_@{cc};

	log(V_@{cc}/vss_@{cc})=rho_v_@{cc}@{cc}*log(V_@{cc}{-1}/vss_@{cc})	
	@#if  strcmp('@{cc}','H')
		+ rho_v_HF*log(V_F{-1}/vss_F)
	@#else
		+ rho_v_FH*log(V_H{-1}/vss_H)
	@#end
	+kappa_v_@{cc}*log(V_HF/v_hf_ss)+sig_v_@{cc}*EV_@{cc};

	log(G_@{cc}/gss_@{cc})=rho_g_@{cc}@{cc}*log(G_@{cc}{-1}/gss_@{cc})	
	@#if  strcmp('@{cc}','H')
		+ rho_g_HF*log(G_F{-1}/gss_F)
	@#else
		+ rho_g_FH*log(G_H{-1}/gss_H)
	@#end
	+sig_g_@{cc}*EG_@{cc};

	R_GC_@{cc} = G_@{cc}/C_@{cc};
@#end

	N_H = P_A*A_F/(
	@#if ncp_shocks
		M_HF{-1}*
	@#end
	V_HF{-1}^(alpha/(1-alpha))*Z_HF{-1})-P_B*B_H;

	N_F = (P_B*B_H*
	@#if ncp_shocks
		M_HF{-1}*
	@#end
	V_HF{-1}^(alpha/(1-alpha))*Z_HF{-1}-P_A*A_F)/P_F;

	RER=P_F;

	TOT=P_B/P_A;

	Y_A=A_H+A_F/(
	@#if ncp_shocks
		M_HF{-1}*
	@#end
	V_HF{-1}^(alpha/(1-alpha))*Z_HF{-1});

	P_H=1;

	% the equation below is lead one step and differs from the one in the note...
	@#if ncp_shocks
		M_HF*
	@#end
	V_HF^(alpha/(1-alpha))*Z_HF*D_H+D_F=0; % updating all the dates from the previous equation

	@#if ncp_shocks
		M_HF=M_H/M_F*M_HF{-1};

		G_LH=L_H/L_H{-1}*M_H{-1};

		R_LFH=L_F/L_H*1/M_HF{-1};
	@#end

	Z_HF=Z_H/Z_F*Z_HF{-1};

	V_HF=V_H/V_F*V_HF{-1};

	G_CH=C_H/C_H{-1}*
		@#if ncp_shocks
			M_H{-1}*
		@#end
		V_H{-1}^(alpha/(1-alpha))*Z_H{-1};

	G_IH=I_H/I_H{-1}*
		@#if ncp_shocks
			M_H{-1}*
		@#end
		V_H{-1}^(1/(1-alpha))*Z_H{-1};

	R_CFH = C_F/C_H*1/(
		@#if ncp_shocks
			M_HF{-1}*
		@#end
		V_HF{-1}^(alpha/(1-alpha))*Z_HF{-1});

	R_IFH =  I_F/I_H*1/(
		@#if ncp_shocks
			M_HF{-1}*
		@#end
		V_HF{-1}^(1/(1-alpha))*Z_HF{-1});

	% Measurement equations
	%-----------------------

	G_CH_HAT=log(G_CH/steady_state(G_CH));
	
	G_IH_HAT=log(G_IH/steady_state(G_IH));
	
	R_GC_H_HAT=log(R_GC_H/steady_state(R_GC_H));
	
	R_CFH_HAT=log(R_CFH/steady_state(R_CFH));
	
	R_IFH_HAT=log(R_IFH/steady_state(R_IFH));
	
	R_GC_F_HAT=log(R_GC_F/steady_state(R_GC_F));
	
	L_H_HAT=log(L_H/steady_state(L_H));
	
	L_F_HAT=log(L_F/steady_state(L_F));
	
	@#if ncp_shocks
		G_LH_HAT=log(G_LH/steady_state(G_LH));
		
		R_LFH_HAT=log(R_LFH/steady_state(R_LFH));
	@#end

