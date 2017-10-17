% Conversion of Dynare file [Dynare_edo.mod] into RISE file [Dynare_edo.rs]
%
% Done 22-Jul-2017 17:34:51.

endogenous

RC	"Marginal product of capital in the slow-growth sector"
RK	"Marginal product of capital in the fast-growth sector"
WC	"Detrended level of wages in the slow-growth sector"
WK	"Detrended level of wages in the fast-growth sector"
YC	"Detrended level of output in the slow-growth sector"
YK	"Detrended level of output in the fast-growth sector"
MCC	"Price markup in the slow-growth sector"
MCK	"Price markup in the fast-growth sector"
KC	"Detrended capital stock in the slow-growth sector"
KK	"Detrended capital stock in the fast-growth sector"
PKB	"Level of replacement price of capital goods, relative to the price of consumer goods"
R	"Federal funds rate"
L	"The marginal utility of the EC concept of consumption"
QK	"Market price of consumer durables"
HC	"Hours in slow-growth sector"
HSC	"Hours in fast-growth sector"
HK	"Potential labor input in the slow-growth sector"
HSK	"Potential labor input in the fast-growth sector"
UHC	"Disutility of labor supply in the slow-growth sector"
UHSC	"Disutility of labor supply potential in the slow-growth sector"
UHK	"Disutility of labor supply in the fast-growth sector"
UHSK	"Disutiltiy of labor potential in the fast-growth sector"
empC	"Employment in the slow-growth sector"
HrC	"Hours per worker in the slow-growth sector"
empK	"Employment in the fast-growth sector"
HrK	"Hours per worker in the fast-growth sector"
empSC	"Potential Employment in the slow-growth sector"
HrSC	"Potential hours per worker in the slow-growth sector"
empSK	"Potential Employment in the fast-growth sector"
HrSK	"Potential hours per worker in the fast-growth sector"
unemp	"Unemployment rate"
EIK	"Detrended level of business investment spending"
EC	"Detrended level of household spending on nondurables and nonhousing services"
INFWC	"Log difference in wages in the slow-growth sector"
INFWK	"Log difference in wages in the fast-growth sector"
INFC	"Consumer-price inflation"
INFK	"Investment-goods inflation"
DIFFNORMGDP	"Log difference of detrended (normalized) real GDP"
NORMINFGDP	"Growth rate of relative price of investment goods"
DIFFREALGDP	"Log difference in real GDP, relative to hours dete"
DIFFREALEC	"Log difference in real household spending on nondurables and nonhousing services (same as FRB/US ECO); relative to hours trend"
DIFFREALEIK	"Log difference in real business spending on investment goods (not BFI, because it includes inventories); relative to hours trend"
DIFFREALW	"Log difference in real household spending on durables (Highly correlated with PIPL-400*d(log(PCNIA)) in FRB/US mnemonics)"
AH	"NFB hours (LHP in FRB/US); log of ratio to trend"
INFGDP	"GDP inflation"
INFCNA	"Overall PCE inflation; PICNIA (/400) in FRB/US"
INFCOR	"Core PCE inflation; PICXFE (/400) in FRB/US"
GAP	"Output gap, Beveridge-Nelson definition"
PFGAP	"Production-function-based output gap"
INFC10	"10 year ahead average core inflation"
ECD	"Detrended level of household spending on durables"
KD	"Detrended stock of consumer durables"
RCD	"Rental rate for consumer durables"
QCD	"Market price of consumer durables"
KCH	"Detrended housing stock"
RCH	"Rental rate for housing"
ECH	"Detrended level of household spending on residential construction"
QCH	"Market price of housing"
LAGKD	"One period lagged detrended stock of consumer durables"
LAGKCH	"One period lagged detrended housing stock"
UK	"Utilization rate of capital in the fast-growth sector"
UC	"Utilization rate of capital in the slow-growth sector"
DIFFREALECH	"Log difference in real household spending on residential construction (same as EH in FRB/US); relative to hours trend"
DIFFREALECD	"Log difference in real business spending on investment goods (not BFI, because it includes inventories); relative to hours trend"
betas	"Aggregate risk premium"
XiL	"Wage-markup shock"
Lpref	"Labor supply preference shock"
EFFK	"Business-investment-specific risk premium"
MUZK	"Growth rate of capital-specific technology shock"
MUZM	"Growth rate of economy-wide technology shock"
HG	"Exogenous spending, normalized"
MUC	"Growth rate of technology in the slow-growth sector"
MUK	"Growth rate of technology in the fast-growth sector"
EFFECD	"Consumer-durables-specific risk premium"
EFFECH	"Housing-specific risk premium"
STAR	"Two-year Treasury term premium"
RL1	"One period federal funds rate lead"
RL2	"Two period federal funds rate lead"
RL3	"Three period federal funds rate lead"
RL4	"Four period federal funds rate lead"
RL5	"Five period federal funds rate lead"
RL6	"Six period federal funds rate lead"
RL7	"Seven period federal funds rate lead"
RT2	"Two-year treasury rate"
DIFFREALGDP_obs	"Observed log difference in real GDP; relative to hours trend"
DIFFREALEC_obs	"Observed log difference in real household spending on nondurables and nonhousing services (same as FRB/US ECO); relative to hours trend"
DIFFREALEIK_obs	"Observed log difference in real business spending on investment goods (not BFI, because it includes inventories); relative to hours trend"
DIFFREALECD_obs	"Observed log difference in real business spending on investment goods (not BFI, because it includes inventories); relative to hours trend"
DIFFREALECH_obs	"Observed log difference in real household spending on residential construction (same as EH in FRB/US); relative to hours trend"
DIFFREALW_obs	"Observed log difference in real household spending on durables (Highly correlated with PIPL-400*d(log(PCNIA)) in FRB/US mnemonics)"
AH_obs	"Observed NFB hours (LHP in FRB/US); log of ratio to trend"
INFCNA_obs	"Observed overall PCE inflation; PICNIA (/400) in FRB/US"
INFCOR_obs	"Observed core PCE inflation; PICXFE (/400) in FRB/US"
INFK_obs	"Observed investment-goods inflation"
R_obs	"Observed federal funds rate"
RT2_obs	"Observed two-year treasury rate"
unemp_obs	"Observed unemployment rate"

exogenous eHG "Exogenous spending shock" eXiL "Wage markup shock" eLpref "The innovation to the labor supply shock" eR "Shock to monetary-policy reaction function"
eMUZK "Shock to investment-specific technology" eMUZM "Shock to economy-wide technology" ePMKC "Price markup shock in slow-growth sector"
ePMKK "Price markup shock in fast-growing sector" eEFFECH "Residential construction risk-premium shock" eEFFECD "Consumer durables risk-premium shock"
eEFFK "Business investment risk-premium shock" eB "Aggregate risk-premium shock" eSTAR "Term-premium shock"


parameters h r_inf r_y r_dy phi_pc phi_H phi_wc phi_ic phi_cd phi_ech gam_pc gam_wc gam_ic gam_icd rho_R rho_B rho_STAR rho_EFFK rho_EFFECD rho_HG rho_EFFECH tp2
ONE MUZMSS MUZKSS r_dinf rpr phi_u rho_MUZK rho_MUZM pbeta delta_ h_cd h_ch delta_cd delta_ch alpha_ theta_c theta_k theta_wc theta_wk g_y a_ks s_AS gam_h
gam_ech s_k s_ecdc eta_cnn eta_cd eta_ch icoef mu_ betarl MUZCSS RCSS RKSS WCSS WKSS YCSS YKSS MCCSS MCKSS KCSS KKSS LSS HCSS HKSS QKSS PKBSS RSS ECSS EIKSS
INFCSS INFKSS INFWCSS INFWKSS MUCSS MUKSS AHSS ECDSS KCDSS QCDSS RCDSS ECHSS KCHSS QCHSS RCHSS UKSS UCSS USS MUKSShabit MUCSShabit INFCNASS INFCORSS INFC10SS
RT2SS beta_0 beta_2 beta_ PYSS AA DD RR eta_cd_eta_cnn eta_ch_eta_cnn Rnr ycbi_ykb hc_hk HSS ycbi ykb YYSS s_k_ecd s_c_ech s_k_eik s_yc IMPHSSS INFGDPSS sig_HG
sig_XiL sig_lpref sig_R sig_MUZK sig_MUZM sig_PMKC sig_PMKK sig_EFFECH sig_EFFECD sig_EFFK sig_B sig_STAR HSKSS HSCSS HrCSS HrKSS A_HC sigman sigmah A_HK xsi_NC
xsi_HrC xsi_NK xsi_HrK rho_XiL rho_lpref empCSS empKSS HrSKSS HrSCSS empSCSS empSKSS UHCSS UHKSS UHSCSS UHSKSS unempSS DIFFREALGDPSS DIFFREALECSS DIFFREALECDSS
DIFFREALECHSS DIFFREALEIKSS DIFFREALWSS RL1SS RL2SS RL3SS RL4SS RL5SS RL6SS RL7SS DIFFREALGDPSS_obs DIFFREALECSS_obs DIFFREALEIKSS_obs DIFFREALECDSS_obs
DIFFREALECHSS_obs DIFFREALWSS_obs INFCNASS_obs INFCORSS_obs INFKSS_obs RSS_obs RT2SS_obs unempSS_obs std_eHG std_eXiL std_eLpref std_eR std_eMUZK std_eMUZM
std_ePMKC std_ePMKK std_eEFFECH std_eEFFECD std_eEFFK std_eB std_eSTAR 

observables DIFFREALGDP_obs DIFFREALEC_obs DIFFREALEIK_obs DIFFREALECD_obs DIFFREALECH_obs DIFFREALW_obs AH_obs INFCNA_obs INFCOR_obs INFK_obs R_obs RT2_obs unemp_obs

model 

RC-MCC*YC/UC/KC(-1)*alpha_*MUK=0;

RK-MCK*YK/UK/KK(-1)*alpha_*MUK=0;

WC-MCC*YC/HC*(1-alpha_)=0;

WK-MCK*YK/HK*(1-alpha_)=0;

YC-(UC*KC(-1)/MUK)^alpha_*(HC)^(1-alpha_)=0;

YK-(UK*KK(-1)/MUK)^alpha_*(HK)^(1-alpha_)=0;

MCC*YC*theta_c-(theta_c-1)*YC-100*phi_pc*(INFC-gam_pc*INFC(-1)-(1-gam_pc)*INFCSS)*INFC*YC+beta_*100*phi_pc*(INFC(+1)-gam_pc*INFC-(1-gam_pc)*INFCSS)*L(+1)/L*INFC(+1)*YC(+1)+100*YCSS*std_ePMKC*ePMKC=0;

MCK*YK*theta_k/PKB-(theta_k-1)*YK-100*phi_pc*(INFK-gam_pc*INFK(-1)-(1-gam_pc)*INFKSS)*INFK*YK+beta_*100*phi_pc*(INFK(+1)-gam_pc*INFK-(1-gam_pc)*INFKSS)*L(+1)/L*YK(+1)*INFK(+1)+100*YKSS*std_ePMKK*ePMKK=0;

QK-beta_*(1/EFFK)*(((1-delta_)*QK(+1)+RC(+1)*UC(+1))*L(+1)/MUK(+1)/L)=0;

QK-beta_*(1/EFFK)*(((1-delta_)*QK(+1)+RK(+1)*UK(+1))*L(+1)/MUK(+1)/L)=0;

L-betas*R/rpr/INFC(+1)/MUC(+1)*L(+1)=0;

log(R/RSS)-rho_R*log(R(-1)/RSS)-(1-rho_R)*(r_inf*log(INFCNA/INFCNASS)+r_dinf*(log(INFCNA)-log(INFCNA(-1)))+r_y*(log(PFGAP)))-std_eR*eR=0;

L-eta_cnn/(EC-h*EC(-1)/MUC)+eta_cnn*beta_*h/(MUC(+1)*EC(+1)-h*EC)=0;

KK-(1-delta_)*KK(-1)/MUK+KC-(1-delta_)*KC(-1)/MUK-1*EIK+mu_*((UK^(1+1/phi_u)-1)/(1+1/phi_u))*KKSS+mu_*((UC^(1+1/phi_u)-1)/(1+1/phi_u))*KCSS=0;

-100+UHC*theta_wc-(theta_wc-1)*WC-100*phi_wc*(INFWC-gam_wc*INFWC(-1)-(1-gam_wc)*INFWCSS)*INFWC*WC+beta_*100*phi_wc*(INFWC(+1)-(gam_wc*INFWC+(1-gam_wc)*INFWCSS))*L(+1)/L*INFWC(+1)*WC(+1)+theta_wc*phi_H/10*(HC/HK-gam_h*HC(-1)/HK(-1)-(1-gam_h)*HCSS/HKSS)+100*XiL=0;

UHSC-WC+phi_H/10*(HSC/HSK-gam_h*HSC(-1)/HSK(-1)-(1-gam_h)*HSCSS/HSKSS);

-100+UHK*theta_wk-(theta_wk-1)*WK-100*phi_wc*(INFWK-gam_wc*INFWK(-1)-(1-gam_wc)*INFWKSS)*INFWK*WK+beta_*100*phi_wc*(INFWK(+1)-(gam_wc*INFWK+(1-gam_wc)*INFWKSS))*L(+1)/L*INFWK(+1)*WK(+1)-theta_wc*phi_H/10*(HC/HK-gam_h*HC(-1)/HK(-1)-(1-gam_h)*HCSS/HKSS)+100*XiL=0;

UHSK-WK-phi_H/10*(HSC/HSK-gam_h*HSC(-1)/HSK(-1)-(1-gam_h)*HSCSS/HSKSS);

UHC*L*Lpref-A_HC*((1+sigman)/(1+sigman/(1+sigmah)))*(HC)^(-1+(1+sigman)/(1+sigman/(1+sigmah)))=0;

UHSC*L*Lpref-A_HC*((1+sigman)/(1+sigman/(1+sigmah)))*(HSC)^(-1+(1+sigman)/(1+sigman/(1+sigmah)))=0;

UHK*L*Lpref-A_HK*((1+sigman)/(1+sigman/(1+sigmah)))*(HK)^(-1+(1+sigman)/(1+sigman/(1+sigmah)))=0;

UHSK*L*Lpref-A_HK*((1+sigman)/(1+sigman/(1+sigmah)))*(HSK)^(-1+(1+sigman)/(1+sigman/(1+sigmah)))=0;

empC-((1+sigmah)/sigmah*xsi_NC/xsi_HrC)^(-1/(1+sigmah+sigman))*HC^(1/(1+sigman/(1+sigmah)))=0;

HrC-((1+sigmah)/sigmah*xsi_NC/xsi_HrC)^(1/(1+sigmah))*empC^(sigman/(1+sigmah))=0;

empK-((1+sigmah)/sigmah*xsi_NK/xsi_HrK)^(-1/(1+sigmah+sigman))*HK^(1/(1+sigman/(1+sigmah)))=0;

HrK-((1+sigmah)/sigmah*xsi_NK/xsi_HrK)^(1/(1+sigmah))*empK^(sigman/(1+sigmah))=0;

empSC-((1+sigmah)/sigmah*xsi_NC/xsi_HrC)^(-1/(1+sigmah+sigman))*HSC^(1/(1+sigman/(1+sigmah)))=0;

HrSC-((1+sigmah)/sigmah*xsi_NC/xsi_HrC)^(1/(1+sigmah))*empSC^(sigman/(1+sigmah))=0;

empSK-((1+sigmah)/sigmah*xsi_NK/xsi_HrK)^(-1/(1+sigmah+sigman))*HSK^(1/(1+sigman/(1+sigmah)))=0;

HrSK-((1+sigmah)/sigmah*xsi_NK/xsi_HrK)^(1/(1+sigmah))*empSK^(sigman/(1+sigmah))=0;

unemp-(empSC+empSK-(empC+empK))/(empSC+empSK)=0;

PKB-(1-100*phi_ic*(EIK-gam_ic*EIK(-1)-(1-gam_ic)*EIKSS)/(KC(-1)+KK(-1))*MUK)*QK-beta_*(1/EFFK)*100*phi_ic*gam_ic*(EIK(+1)-(gam_ic*EIK+(1-gam_ic)*EIKSS))/(KC+KK)*MUK(+1)*QK(+1)*L(+1)/L=0;

YC-EC-ECH-0.2*YCSS*HG=0;

log(INFWC)-log(WC)+log(WC(-1))-log(MUC)-log(INFC)=0;

log(INFWK)-log(WK)+log(WK(-1))-log(MUC)-log(INFC)=0;

log(INFK)-log(INFC)-log(PKB)+log(PKB(-1))+log(MUK)-log(MUC)=0;

YK-EIK-ECD-0.2*YKSS*HG=0;

log(DIFFNORMGDP)-(1-s_k)*(log(YC)-log(YC(-1)))-s_k*(log(YK)-log(YK(-1)))=0;

log(NORMINFGDP)-s_k*(log(PKB)-log(PKB(-1)))=0;

log(DIFFREALGDP)-log(DIFFNORMGDP)-(1-s_k)*log(MUC)-s_k*log(MUK)=0;

log(DIFFREALEC)-log(EC)+log(EC(-1))-log(MUC)=0;

log(DIFFREALEIK)-log(EIK)+log(EIK(-1))-log(MUK)=0;

log(DIFFREALW)-HCSS/AHSS*(log(INFWC))-HKSS/AHSS*(log(INFWK))+log(INFC)=0;

AH-HC-HK=0;

log(INFGDP)-log(INFC)-log(YC*MUC/YC(-1))+log(DIFFREALGDP)-log((1+PKB*YK/YC)/(1+PKB(-1)*YK(-1)/YC(-1)))=0;

log(INFCNA)-(1-s_ecdc)*log(INFC)-s_ecdc*log(INFK)=0;

log(INFCOR)-(1-s_ecdc)*log(INFC)-s_ecdc*log(INFK)=0;

log(GAP)-(1-s_k)*log(YC/YCSS)-s_k*log(YK/YKSS)=0;

log(PFGAP)-(1-alpha_)*((1-s_k)*log(HC/HCSS)+s_k*log(HK/HKSS))-alpha_*((1-s_k)*log(UC/USS)+s_k*log(UK/USS))=0;

log(INFC10)-betarl*log(INFC10(+1))-(1-betarl)*log(INFCOR)=0;

KD-(1-delta_cd)*KD(-1)/MUK-ECD=0;

L*RCD-eta_cd/(KD(-1)/MUK-h_cd*LAGKD(-1)/(MUK(-1)*MUK))+beta_*eta_cd*h_cd/(KD-h_cd*KD(-1)/MUK)=0;

QCD-beta_*(1/EFFECD)*L(+1)/L/MUK(+1)*(RCD(+1)+(1-delta_cd)*QCD(+1))=0;

PKB-QCD*(1-100*phi_cd*(ECD-gam_icd*ECD(-1)-(1-gam_icd)*ECDSS)/KD(-1)*MUK) - beta_*(1/EFFECD)*100*gam_icd*phi_cd*(ECD(+1)-(gam_icd*ECD+(1-gam_icd)*ECDSS))/KD*QCD(+1)*L(+1)/L*MUK(+1)=0;

L*RCH-eta_ch/(KCH(-1)/MUC-h_ch*LAGKCH(-1)/(MUC*MUC(-1)))+beta_*eta_ch*h_ch/(KCH-h_ch*KCH(-1)/MUC)=0;

QCH-beta_*(1/EFFECH)*L(+1)/L/MUC(+1)*(RCH(+1)+(1-delta_ch)*QCH(+1))=0;

1*ECH+(1-delta_ch)*KCH(-1)/MUC-KCH=0;

1-QCH*(1-100*phi_ech*(ECH-gam_ech*ECH(-1)-(1-gam_ech)*ECHSS)/KCH(-1)*MUC) - beta_*(1/EFFECH)*100*gam_ech*phi_ech*(ECH(+1)-(gam_ech*ECH+(1-gam_ech)*ECHSS))/KCH*QCH(+1)*L(+1)/L*MUC(+1)=0;

log(KD(-1))-log(LAGKD)=0;

log(KCH(-1))-log(LAGKCH)=0;

RK-QK*mu_*UK^(1/phi_u)=0;

RC-QK*mu_*UC^(1/phi_u)=0;

log(DIFFREALECH)-log(MUC)-log(ECH)+log(ECH(-1))=0;

log(DIFFREALECD)-log(MUK)-log(ECD)+log(ECD(-1))=0;

log(betas/beta_)-rho_B*log(betas(-1)/beta_)-std_eB*eB=0;

log(XiL)-rho_XiL*log(XiL(-1))-std_eXiL*eXiL=0;

log(Lpref)-rho_lpref*log(Lpref(-1))-std_eLpref*eLpref=0;

log(EFFK)-rho_EFFK*log(EFFK(-1))-std_eEFFK*eEFFK=0;

log(MUZK/MUZKSS)-std_eMUZK*eMUZK=0;

log(MUZM/MUZMSS)-std_eMUZM*eMUZM=0;

log(HG)-rho_HG*log(HG(-1))-std_eHG*eHG=0;

log(MUC)-log(MUZM)-alpha_*log(MUZK)=0;

log(MUK)-log(MUZM)-log(MUZK)=0;

log(EFFECD)-rho_EFFECD*log(EFFECD(-1))-std_eEFFECD*eEFFECD=0;

log(EFFECH)-rho_EFFECH*log(EFFECH(-1))-std_eEFFECH*eEFFECH=0;

log(STAR)-rho_STAR*log(STAR(-1))-std_eSTAR*eSTAR=0;

log(RL1) - log(R(+1))=0;

log(RL2) - log(RL1(+1))=0;

log(RL3) - log(RL2(+1))=0;

log(RL4) - log(RL3(+1))=0;

log(RL5) - log(RL4(+1))=0;

log(RL6) - log(RL5(+1))=0;

log(RL7) - log(RL6(+1))=0;

log(RT2) - tp2 - 0.125*(log(R) + log(RL1) + log(RL2) + log(RL3) + log(RL4) + log(RL5) + log(RL6) + log(RL7)) - log(STAR)=0;

log(DIFFREALGDP_obs/DIFFREALGDPSS_obs) = log(DIFFREALGDP/DIFFREALGDPSS);

log(DIFFREALEC_obs/DIFFREALECSS_obs) = log(DIFFREALEC/DIFFREALECSS);

log(DIFFREALEIK_obs/DIFFREALEIKSS_obs) = log(DIFFREALEIK/DIFFREALEIKSS);

log(DIFFREALECD_obs/DIFFREALECDSS_obs) = log(DIFFREALECD/DIFFREALECDSS);

log(DIFFREALECH_obs/DIFFREALECHSS_obs) = log(DIFFREALECH/DIFFREALECHSS);

log(DIFFREALW_obs/DIFFREALWSS_obs) = log(DIFFREALW/DIFFREALWSS);

log(AH_obs) = log(AH/AHSS);

log(INFCNA_obs/INFCNASS_obs) = log(INFCNA/INFCNASS);

log(INFCOR_obs/INFCORSS_obs) = log(INFCOR/INFCORSS);

log(INFK_obs/INFKSS_obs) = log(INFK/INFKSS);

log(R_obs/RSS_obs) = log(R/RSS);

log(RT2_obs/RT2SS_obs) = log(RT2/RT2SS);

log(unemp_obs/unempSS_obs) = log(unemp/unempSS);



