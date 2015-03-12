% model1.rs line 9
endogenous
C_H, "Consumption(H)" I_H, "Investment(H)" G_H,  "Gov. Spending process(H)" N_H, "Net Exports(H)"
CTILDE_H, "Total consumption(H)" ITILDE_H, "Total investment(H)" K_H, "Capital(H)"
D_H, "Bonds(H)" LAMBDA_H "Lagrange mult. on budget(H)" XI_H, "Lagrange mult. on capital accum.(H)"
TAU_H, "Lump-sum taxes(H)" A_H, "Demand of home interm. goods(H)" 
B_H, "Demand of Foreign interm. goods(H)" P_H, "Price of cons. goods(H)" 
X_H, "Price of Inv. goods(H)" W_H, "Nom. wages(H)" 
Q_H, "Nom. rental rate of capital(H)" M_H, "Preference shock process(H)" 
Z_H, "Neutral tech. shock process(H)", V_H "Inv.-specific tech. shock process(H)"
L_H, "Hours worked(H)" R_GC_H "Gov. spending-Consumption ratio(H)" 
C_F, "Consumption(F)" I_F, "Investment(F)" G_F,  "Gov. Spending process(F)" N_F, "Net Exports(F)"
CTILDE_F, "Total consumption(F)" ITILDE_F, "Total investment(F)" K_F, "Capital(F)"
D_F, "Bonds(F)" LAMBDA_F "Lagrange mult. on budget(F)" XI_F, "Lagrange mult. on capital accum.(F)"
TAU_F, "Lump-sum taxes(F)" A_F, "Demand of home interm. goods(F)" 
B_F, "Demand of Foreign interm. goods(F)" P_F, "Price of cons. goods(F)" 
X_F, "Price of Inv. goods(F)" W_F, "Nom. wages(F)" 
Q_F, "Nom. rental rate of capital(F)" M_F, "Preference shock process(F)" 
Z_F, "Neutral tech. shock process(F)", V_F "Inv.-specific tech. shock process(F)"
L_F, "Hours worked(F)" R_GC_F "Gov. spending-Consumption ratio(F)" 
Y_A, "Home intermediate goods" Y_B "Foreign intermediate goods", R "Interest rate",
P_A, "real price of home interm. goods" P_B, "real price of foreign interm. goods" RER, "real exch. rate (pF)"
TOT, "Terms of trade(pB/pA)" Z_HF, "Neutral tech. ratio(H/F)" V_HF, "Inv.-specific tech. ratio(H/F)"
G_CH, "Consumption growth"
G_IH, "Inv. growth" R_CFH, "Foreign-Home consumption ratio" R_IFH "Foreign/Home Invest. ratio"
G_CH_HAT, G_IH_HAT, R_GC_H_HAT, L_H_HAT, R_CFH_HAT, R_IFH_HAT, R_GC_F_HAT, L_F_HAT

observables G_CH_HAT, G_IH_HAT, R_GC_H_HAT, L_H_HAT, R_CFH_HAT, R_IFH_HAT, R_GC_F_HAT, L_F_HAT

exogenous
	EM_H "Preference(H)" EV_H "Inv.-specific tech(H)" EZ_H "Neutral tech(H)" EG_H "Gov. spending(H)"
	EM_F "Preference(F)" EV_F "Inv.-specific tech(F)" EZ_F "Neutral tech(F)" EG_F "Gov. spending(F)"
	
parameters
alpha, phi_d philtr phiktr mu gam delta theta omega z_hf_ss v_hf_ss beta kappazhtr
kappavhtr
eta_H,
rho_m_HH, rho_z_HH, rho_v_HH, rho_g_HH, 
	 rho_m_HF,rho_z_HF, rho_v_HF, rho_g_HF, 
mss_H, zss_H, vss_H, gss_H
sig_g_H	sig_v_H	sig_z_H sig_m_H
eta_F,
rho_m_FF, rho_z_FF, rho_v_FF, rho_g_FF, 
	rho_m_FH, rho_z_FH, rho_v_FH, rho_g_FH
	kappa_z_F, kappa_v_F
mss_F, zss_F, vss_F, gss_F
sig_g_F	sig_v_F	sig_z_F sig_m_F
	
model
	# phi_l = 100*abs(philtr);
	# phi_k = 100*abs(phiktr);
	# kappa_z_H = -abs(kappazhtr);
	# kappa_v_H = -abs(kappavhtr);
	
	W_H*L_H + Q_H*K_H{-1} + D_H{-1} = C_H + X_H*I_H + TAU_H + V_H^(alpha/(1-alpha))*Z_H*D_H/(R)
		+phi_d/2*V_H^(alpha/(1-alpha))*Z_H*D_H^2+phi_l/2*(L_H/L_H{-1}-1)^2*L_H{-1};	% (1)
	
	(1-delta)*K_H{-1}+I_H-phi_k/2*(I_H/K_H{-1}-eta_H)^2*K_H{-1}=V_H^(1/(1-alpha))*Z_H*K_H;	% (2)
	
	mu*(C_H^mu*(1-L_H/M_H)^(1-mu))^(1-gam)=LAMBDA_H*C_H;	% (3)
	
	(1-mu)*(C_H^mu*(1-L_H/M_H)^(1-mu))^(1-gam)/(M_H*(1-L_H/M_H))= LAMBDA_H*(W_H-phi_l*(L_H/L_H{-1}-1))+
		beta*phi_l*LAMBDA_H{+1}*(V_H^(alpha/(1-alpha))*Z_H)^(mu*(1-gam))*((L_H{+1}/L_H-1)*L_H{+1}/L_H
		-1/2*(L_H{+1}/L_H-1)^2);	% (4)

	LAMBDA_H*X_H = XI_H*(1-phi_k*(I_H/K_H{-1}-eta_H));	% (5)
	
	(V_H^(alpha/(1-alpha))*Z_H)^(1-mu*(1-gam))*V_H*XI_H = beta*LAMBDA_H{+1}*Q_H{+1} +
		beta*XI_H{+1}*(1-delta+phi_k*(I_H{+1}/K_H-eta_H)*I_H{+1}/K_H-phi_k/2*(I_H{+1}/K_H-eta_H)^2);	% (6)
		
	(V_H^(alpha/(1-alpha))*Z_H)^(1-mu*(1-gam))*LAMBDA_H*(1/(R)+phi_d*D_H)=beta*LAMBDA_H{+1};	% (7)
	
	W_F*L_F+Q_F*K_F{-1}+D_F{-1}/P_F = C_F+X_F*I_F+TAU_F+V_F^(alpha/(1-alpha))*Z_F*D_F/(
		P_F*R)+phi_d/2*V_F^(alpha/(1-alpha))*Z_F*D_F^2+phi_l/2*(L_F/L_F{-1}-1)^2*L_F{-1};	% (8)
		
	(1-delta)*K_F{-1}+I_F-phi_k/2*(I_F/K_F{-1}-eta_F)^2*K_F{-1}=V_F^(1/(1-alpha))*Z_F*K_F;	% (9)
	
	mu*(C_F^mu*(1-L_F/M_F)^(1-mu))^(1-gam)=LAMBDA_F*C_F;	% (10)
	
	(1-mu)*(C_F^mu*(1-L_F/M_F)^(1-mu))^(1-gam)/(M_F*(1-L_F/M_F))= LAMBDA_F*(W_F-phi_l*(L_F/L_F{-1}-1))+
		beta*phi_l*LAMBDA_F{+1}*(V_F^(alpha/(1-alpha))*Z_F)^(mu*(1-gam))*((L_F{+1}/L_F-1)*L_F{+1}/L_F
		-1/2*(L_F{+1}/L_F-1)^2);	% (11)
		
	LAMBDA_F*X_F = XI_F*(1-phi_k*(I_F/K_F{-1}-eta_F));	% (12)
	
	(V_F^(alpha/(1-alpha))*Z_F)^(1-mu*(1-gam))*V_F*XI_F = beta*LAMBDA_F{+1}*Q_F{+1} +
		beta*XI_F{+1}*(1-delta+phi_k*(I_F{+1}/K_F-eta_F)*I_F{+1}/K_F
		-phi_k/2*(I_F{+1}/K_F-eta_F)^2);	% (13)
		
	(V_F^(alpha/(1-alpha))*Z_F)^(1-mu*(1-gam))*LAMBDA_F*(1/(P_F*R)+phi_d*D_F)=beta*LAMBDA_F{+1}/P_F{+1};	% (14)
	
	K_H{-1}^alpha*(Z_H*L_H)^(1-alpha)=Y_A;	% (15)
	
	alpha*P_A*Y_A = Q_H*K_H{-1};	% (16)
	
	(1-alpha)*P_A*Y_A=W_H*L_H;	% (17)

	K_F{-1}^alpha*(Z_F*L_F)^(1-alpha)=Y_B;	% (18)
	
	alpha*P_B*Y_B=P_F*Q_F*K_F{-1};	% (19)
		
	(1-alpha)*P_B*Y_B=P_F*W_F*L_F;	% (20)
		
	((1-omega)^(1/theta)*A_H^((theta-1)/theta)+omega^(1/theta)*B_H^((theta-1)/theta))^(theta/(theta-1))=
		CTILDE_H+1/V_H*ITILDE_H;	% (21)
	
	X_H = 1/V_H;	% (22)
	
	P_A= ((1-omega)^(1/theta)*A_H^((theta-1)/theta)+omega^(1/theta)*B_H^((theta-1)/theta))^(1/(theta-1))*(1-omega)
		^(1/theta)*A_H^(-1/theta);	% (23)
	
	P_B= ((1-omega)^(1/theta)*A_H^((theta-1)/theta)+omega^(1/theta)*B_H^((theta-1)/theta))^(1/(theta-1))*
		omega^(1/theta)*B_H^(-1/theta);	% (24)
		
	(omega^(1/theta)*A_F^((theta-1)/theta)+(1-omega)^(1/theta)*B_F^((theta-1)/theta))^(theta/(theta-1))=
		CTILDE_F+1/V_F*ITILDE_F;	% (25)
		
	X_F = 1/V_F;	% (26)
	
	P_A/P_F= (omega^(1/theta)*A_F^((theta-1)/theta)+(1-omega)^(1/theta)*B_F^((theta-1)/theta))^(1/(theta-1))*
		omega^(1/theta)*A_F^(-1/theta);	% (27)
	
	P_B/P_F= (omega^(1/theta)*A_F^((theta-1)/theta)+(1-omega)^(1/theta)*B_F^((theta-1)/theta))^(1/(theta-1))*
		(1-omega)^(1/theta)*B_F^(-1/theta);	% (28)
		
	TAU_H = G_H;	% (29)
	
	TAU_F = G_F;	% (30)
	
	N_H = P_A*A_F/(V_HF{-1}^(alpha/(1-alpha))*Z_HF{-1})-P_B*B_H;	% (31)
	
	N_F = (P_B*B_H*V_HF{-1}^(alpha/(1-alpha))*Z_HF{-1}-P_A*A_F)/P_F;	% (32)
	
	RER=P_F;	% (33)
	
	TOT=P_B/P_A;	% (34)
	
	Y_A=A_H+A_F/(V_HF{-1}^(alpha/(1-alpha))*Z_HF{-1});	% (35)
	
	P_H=1;	% (36)
	
	CTILDE_H = C_H+phi_d/2*V_H^(alpha/(1-alpha))*Z_H*D_H^2+phi_l/2*(L_H/L_H{-1}-1)^2*L_H{-1}+G_H;	% (37)
	
	ITILDE_H = I_H;	% (38)
	
	CTILDE_F = C_F+phi_d/2*V_F^(alpha/(1-alpha))*Z_F*D_F^2+phi_l/2*(L_F/L_F{-1}-1)^2*L_F{-1}+G_F;	% (39)
	
	ITILDE_F = I_F;	% (40)
	
	V_HF^(alpha/(1-alpha))*Z_HF*D_H+D_F=0;	% (41) %<--V_HF{-1}^(alpha/(1-alpha))*Z_HF{-1}*D_H+D_F=0;
	
	log(M_H/mss_H)=rho_m_HH*log(M_H{-1}/mss_H)+ rho_m_HF*log(M_F{-1}/mss_F)+sig_m_H*EM_H;	% (42)
	
	log(M_F/mss_F)=rho_m_FF*log(M_F{-1}/mss_F) + rho_m_FH*log(M_H{-1}/mss_H)+sig_m_F*EM_F;	% (43)
	
	log(Z_H/zss_H)=rho_z_HH*log(Z_H{-1}/zss_H)+ rho_z_HF*log(Z_F{-1}/zss_F)+kappa_z_H*log(Z_HF/z_hf_ss)+sig_z_H*EZ_H;	% (44)
	
	log(Z_F/zss_F)=rho_z_FF*log(Z_F{-1}/zss_F) + rho_z_FH*log(Z_H{-1}/zss_H) +kappa_z_F*log(Z_HF/z_hf_ss)+sig_z_F*EZ_F;	% (45)
	
	log(V_H/vss_H)=rho_v_HH*log(V_H{-1}/vss_H)+ rho_v_HF*log(V_F{-1}/vss_F)+kappa_v_H*log(V_HF/v_hf_ss)+sig_v_H*EV_H;	% (46)
	
	log(V_F/vss_F)=rho_v_FF*log(V_F{-1}/vss_F) + rho_v_FH*log(V_H{-1}/vss_H) +kappa_v_F*log(V_HF/v_hf_ss)+sig_v_F*EV_F;	% (47)
	
	log(G_H/gss_H)=rho_g_HH*log(G_H{-1}/gss_H)+ rho_g_HF*log(G_F{-1}/gss_F)+sig_g_H*EG_H;	% (48)
	
	log(G_F/gss_F)=rho_g_FF*log(G_F{-1}/gss_F) + rho_g_FH*log(G_H{-1}/gss_H) +sig_g_F*EG_F;	% (49)
	
	Z_HF=Z_H/Z_F*Z_HF{-1};	% (50)
	
	V_HF=V_H/V_F*V_HF{-1};	% (51)
	
	G_CH=C_H/C_H{-1}*V_H{-1}^(alpha/(1-alpha))*Z_H{-1};	% (52)
	
	G_IH=I_H/I_H{-1}*V_H{-1}^(1/(1-alpha))*Z_H{-1}	% (53)
	
	R_GC_H = G_H/C_H;	% (54)
	
	R_CFH = C_F/C_H*1/(V_HF{-1}^(alpha/(1-alpha))*Z_HF{-1});	% (55)
	
	R_IFH =  I_F/I_H*1/(V_HF{-1}^(1/(1-alpha))*Z_HF{-1});	% (56)
	
	R_GC_F = G_F/C_F;	% (57)
	
	G_CH_HAT=log(G_CH/steady_state(G_CH));
	
	G_IH_HAT=log(G_IH/steady_state(G_IH));
	
	R_GC_H_HAT=log(R_GC_H/steady_state(R_GC_H));
	
	L_H_HAT=log(L_H/steady_state(L_H));
	
	R_CFH_HAT=log(R_CFH/steady_state(R_CFH));
	
	R_IFH_HAT=log(R_IFH/steady_state(R_IFH));
	
	R_GC_F_HAT=log(R_GC_F/steady_state(R_GC_F));
	
	L_F_HAT=log(L_F/steady_state(L_F));
