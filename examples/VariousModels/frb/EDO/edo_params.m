function [p,priors]=edo_params()
% Conversion of Dynare file [Dynare_edo.mod] into RISE file [edo.rs]
%
% Done 22-Jul-2017 19:19:44.
% 
% Remarks: 
% - The parameters and shocks variances or 
%   standard deviations found in the dynare file are assigned. 
% - Shock standard deviations not assigned in the 
%   dynare file get a value of 0 following dynare's convention
% - RISE will set all other parameters without a value to nan

p=struct();

p.h=0.715162417869797;
p.r_inf=1.46344163969035;
p.r_y=0.263123294207851;
p.phi_pc=3.54471453295450;
p.phi_H=3.22894079106560;
p.phi_wc=5.49395755514723;
p.phi_ic=0.253308786976374;
p.phi_cd=0.470089385005009;
p.phi_ech=9.13986886546163;
p.gam_pc=0.314488926051065;
p.gam_wc=-0.230018833252054;
p.sigman=39.4075260618789;
p.sigmah=21.8859803402692;
p.rho_R=0.833200065745674;
p.rho_XiL=0.263567746111198;
p.rho_lpref=0.979092048897712;
p.rho_B=0.895267027146152;
p.rho_STAR=0.909187927454138;
p.rho_EFFK=0.937829274540004;
p.rho_EFFECD=-0.240286975088701;
p.rho_HG=0.582395471123139;
p.rho_EFFECH=0.877235725078934;
p.tp2=0.000307314910763576;
p.sig_HG=0.579315931803017;
p.sig_XiL=2.49313873916751;
p.sig_lpref=5.66476748114241;
p.sig_R=0.124100461010359;
p.sig_MUZK=0.936167718269030;
p.sig_MUZM=0.597390920898135;
p.sig_PMKC=0.451830653200989;
p.sig_PMKK=0.685376191952156;
p.sig_EFFECH=0.514704527091087;
p.sig_EFFECD=9.11199585973990;
p.sig_EFFK=0.402779878811407;
p.sig_B=0.295232712196573;
p.sig_STAR=0.104877885500673;
p.r_dy=0;
p.ONE=1;
p.MUZKSS=1.009250;
p.MUZMSS=1.001000;
p.gam_ic=1.0;
p.gam_icd=1.0;
p.r_dinf=0;
p.rpr=0.965;
p.phi_u=1;
p.rho_MUZK=0;
p.rho_MUZM=0;
p.pbeta=0.99862;
p.delta_=0.03;
p.h_cd=0.0;
p.h_ch=0.0;
p.delta_cd=0.055;
p.delta_ch=0.0035;
p.alpha_=0.26;
p.theta_c=7;
p.theta_k=7;
p.unempSS=.06;
p.g_y=0.0;
p.a_ks=0.2;
p.s_AS=0.2;
p.gam_h=1;
p.gam_ech=1;
p.icoef=3;
p.betarl=.958;
p.std_eHG=p.sig_HG;
p.std_eXiL=p.sig_XiL;
p.std_eLpref=p.sig_lpref;
p.std_eR=p.sig_R;
p.std_eMUZK=p.sig_MUZK;
p.std_eMUZM=p.sig_MUZM;
p.std_ePMKC=p.sig_PMKC;
p.std_ePMKK=p.sig_PMKK;
p.std_eEFFECH=p.sig_EFFECH;
p.std_eEFFECD=p.sig_EFFECD;
p.std_eEFFK=p.sig_EFFK;
p.std_eB=p.sig_B;
p.std_eSTAR=p.sig_STAR;




priors=struct();

priors.h={.673, -1, 1};

priors.r_inf={1.461, 1.5000, 0.0625000, 'normal', -999, 999};

priors.r_y={0.214, 0.125, 0.125000, 'normal', -999, 999};

priors.phi_pc={3.126, 4.0000, 4.0000^.5, 'gamma', 0, 999};

priors.phi_H={4.064, 4.0000, 4.0000^.5, 'gamma', 0, 999};

priors.phi_wc={5.119, 4.0000, 4.0000^.5, 'gamma', 0, 999};

priors.phi_ic={.325, 4.0000, 4.0000^.5, 'gamma', 0, 999};

priors.phi_cd={.651, 4.0000, 4.0000^.5, 'gamma', 0, 999};

priors.phi_ech={10.948, 4.0000, 4.0000^.5, 'gamma', 0, 999};

priors.gam_pc={0.386, 0.000, 0.250, 'normal', -999, 999};

priors.gam_wc={0.213, 0.000, 0.250, 'normal', -999, 999};

priors.sigman={1.25, 1.25, 12.5^.5, 'gamma', 0, 999};

priors.sigmah={10, 10, 100^.5, 'gamma', 0, 999};

priors.rho_R={0.654, 0.5, 0.25, 'normal', -1, 1};

priors.rho_XiL={0.654, 0.5, 0.25, 'normal', -1, 1};

priors.rho_lpref={0.954, 0.5, 0.25, 'normal', -1, 1};

priors.rho_B={0.825, 0, 0.5, 'normal', -1, 1};

priors.rho_STAR={0.825, 0, 0.5, 'normal', -1, 1};

priors.rho_EFFK={0.850, 0, 0.5, 'normal', -1, 1};

priors.rho_EFFECD={.230, 0, 0.5, 'normal', -1, 1};

priors.rho_HG={0.596, 0.5, 0.015^.5, 'beta'};

priors.rho_EFFECH={0.844, 0, 0.5, 'normal', -1, 1};

priors.tp2={0.001, 0.0, 0.0005, 'normal', -999, 999};

priors.std_eHG={.745, 1.772454, 4, 'inv_gamma', 0.0001, 999};

priors.std_eXiL={3.621, 1.772454, 4, 'inv_gamma', 0.0001, 999};

priors.std_eLpref={1.621, 1.772454, 4, 'inv_gamma', 0.0001, 999};

priors.std_eR={0.165, 0.354491, 4, 'inv_gamma', 0.0001, 999};

priors.std_eMUZK={.834, 0.443113, 4, 'inv_gamma', 0.0001, 999};

priors.std_eMUZM={.484, 0.443113, 4, 'inv_gamma', 0.0001, 999};

priors.std_ePMKC={.391, 0.354491, 4, 'inv_gamma', 0.0001, 999};

priors.std_ePMKK={.552, 0.354491, 4, 'inv_gamma', 0.0001, 999};

priors.std_eEFFECH={.526, 1.772454, 4, 'inv_gamma', 0.0001, 999};

priors.std_eEFFECD={13.349, 1.772454, 4, 'inv_gamma', 0.0001, 999};

priors.std_eEFFK={.499, 1.772454, 4, 'inv_gamma', 0.0001, 999};

priors.std_eB={0.5, 1.772454, 4, 'inv_gamma', 0.0001, 999};

priors.std_eSTAR={0.05, 0.354491, 4, 'inv_gamma', 0.0001, 999};

