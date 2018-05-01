function [y,newp,retcode]=ssfile(obj,y,p,d,id) %#ok<INUSL>

% risessfile --  computes the steady state of ... analytically
%
% ::
%
%
%   [y,newp,retcode]=risessfile(obj,y,p,d,id)
%
% Args:
%
%    - **obj** [rise|dsge]: model object (not always needed)
%
%    - **y** [vector]: endo_nbr x 1 vector of initial steady state
%
%    - **p** [struct]: parameter structure
%
%    - **d** [struct]: definitions
%
%    - **id** [vector]: location of the variables to calculate
%
% Returns:
%    :
%
%    - **y** []: endo_nbr x 1 vector of updated steady state
%
%    - **newp** [struct]: structure containing updated parameters if any
%
%    - **retcode** [0|number]: return 0 if there are no problems, else return
%      any number different from 0
%
% Note:
%
%    - this is new approach has three main advantages relative to the previous
%      one:
%      - The file is valid whether we have many regimes or not
%      - The user does not need to know what regime is being computed
%      - It is in sync with the steady state model
%
% Example:
%
%    See also:

% Conversion of Dynare file [Dynare_edo_steadystate.m] into RISE file [risessfile.m]
%
%  Done 23-Jul-2017 00:31:15.

retcode=0;

if nargin==1
    % list of endogenous variables to be calculated
    %----------------------------------------------
    % y=get(obj,'endo_list');
    y={'RC','RK','WC','WK','YC','YK','MCC','MCK','KC','KK','PKB','R','L','QK',...
        'HC','HSC','HK','HSK','UHC','UHSC','UHK','UHSK','empC','HrC',...
        'empK','HrK','empSC','HrSC','empSK','HrSK','unemp','EIK','EC','INFWC',...
        'INFWK','INFC','INFK','DIFFNORMGDP','NORMINFGDP','DIFFREALGDP',...
        'DIFFREALEC','DIFFREALEIK','DIFFREALW','AH','INFGDP','INFCNA','INFCOR',...
        'GAP','PFGAP','INFC10','ECD','KD','RCD','QCD','KCH','RCH','ECH','QCH',...
        'LAGKD','LAGKCH','UK','UC','DIFFREALECH','DIFFREALECD','betas','XiL',...
        'Lpref','EFFK','MUZK','MUZM','HG','MUC','MUK','EFFECD','EFFECH','STAR',...
        'RL1','RL2','RL3','RL4','RL5','RL6','RL7','RT2','DIFFREALGDP_obs',...
        'DIFFREALEC_obs','DIFFREALEIK_obs','DIFFREALECD_obs','DIFFREALECH_obs',...
        'DIFFREALW_obs','AH_obs','INFCNA_obs','INFCOR_obs','INFK_obs','R_obs',...
        'RT2_obs','unemp_obs'};
    % list of parameters to be computed during steady state calculation
    %-------------------------------------------------------------------
    newp={'AA','AHSS','A_HC','A_HK','DD','DIFFREALECDSS','DIFFREALECDSS_obs',...
        'DIFFREALECHSS','DIFFREALECHSS_obs','DIFFREALECSS','DIFFREALECSS_obs',...
        'DIFFREALEIKSS','DIFFREALEIKSS_obs','DIFFREALGDPSS','DIFFREALGDPSS_obs',...
        'DIFFREALWSS','DIFFREALWSS_obs','ECDSS','ECHSS','ECSS','EIKSS','HCSS',...
        'HKSS','HSCSS','HSKSS','HSS','HrCSS','HrKSS','HrSCSS','HrSKSS',...
        'IMPHSSS','INFC10SS','INFCNASS','INFCNASS_obs','INFCORSS','INFCORSS_obs',...
        'INFCSS','INFGDPSS','INFKSS','INFKSS_obs','INFWCSS','INFWKSS','KCDSS',...
        'KCHSS','KCSS','KKSS','LSS','MCCSS','MCKSS','MUCSS','MUCSShabit','MUKSS',...
        'MUKSShabit','MUZCSS','ONE','PKBSS','PYSS','QCDSS','QCHSS','QKSS',...
        'RCDSS','RCHSS','RCSS','RKSS','RL1SS','RL2SS','RL3SS','RL4SS','RL5SS',...
        'RL6SS','RL7SS','RR','RSS','RSS_obs','RT2SS','RT2SS_obs','Rnr','UCSS',...
        'UHCSS','UHKSS','UHSCSS','UHSKSS','UKSS','USS','WCSS','WKSS','YCSS',...
        'YKSS','YYSS','beta_','beta_0','beta_2','empCSS','empKSS','empSCSS',...
        'empSKSS','eta_cd','eta_cd_eta_cnn','eta_ch','eta_ch_eta_cnn','eta_cnn',...
        'hc_hk','mu_','s_c_ech','s_ecdc','s_k','s_k_ecd','s_k_eik','s_yc',...
        'theta_wc','theta_wk','unempSS','unempSS_obs','xsi_HrC','xsi_HrK',...
        'xsi_NC','xsi_NK','ycbi','ycbi_ykb','ykb'};
else
    
    newp=struct();
    
    %start_steady_state;
    
    newp.beta_0 = p.pbeta;
    newp.beta_2 = p.pbeta*p.rpr; % s.s. funds rate premium
    newp.beta_ = newp.beta_2;
    newp.MUZCSS=1;
    newp.ONE=1;
    newp.USS=1;
    newp.MUKSS=p.MUZKSS*p.MUZMSS;
    newp.MUCSS=p.MUZKSS^p.alpha_*p.MUZMSS;
    newp.MUKSShabit=newp.MUKSS;
    newp.MUCSShabit=newp.MUCSS;
    newp.PKBSS=p.theta_k/(p.theta_k-1)*(p.theta_c-1)/p.theta_c;
    newp.PYSS=1;
    newp.MCCSS=(p.theta_c-1)/p.theta_c;
    newp.MCKSS=(p.theta_k-1)/p.theta_k;
    newp.RKSS=newp.MUKSS/newp.beta_2-(1-p.delta_);
    newp.RCSS=newp.MUKSS/newp.beta_2-(1-p.delta_);
    newp.RCHSS=newp.MUCSS/newp.beta_2-(1-p.delta_ch); % Housing sector
    newp.RCDSS=newp.MUKSS/newp.beta_2-(1-p.delta_cd); % Durable sector
    newp.USS=1;
    newp.mu_=newp.RCSS;
    newp.AA=p.alpha_/newp.RKSS*newp.MCKSS;
    newp.DD = 0.135;
    newp.RR = 0.075;
    newp.eta_cnn=1;
    newp.eta_cd_eta_cnn=newp.DD/((newp.MUKSShabit-newp.beta_2*p.h_cd)/(1-newp.beta_2*p.h/newp.MUCSShabit)*(1-p.h/newp.MUCSShabit)/(1-p.h_cd/newp.MUKSShabit)*(1-(1-p.delta_cd)/newp.MUKSS)/newp.RCDSS);
    newp.eta_ch_eta_cnn=newp.RR/((newp.MUCSShabit-newp.beta_2*p.h_ch)/(1-newp.beta_2*p.h/newp.MUCSShabit)*(1-p.h/newp.MUCSShabit)/(1-p.h_ch/newp.MUCSShabit)*(1-(1-p.delta_ch)/newp.MUCSS)/newp.RCHSS);
    newp.eta_ch=newp.eta_ch_eta_cnn;
    newp.eta_cd=newp.eta_cd_eta_cnn;
    newp.DD=newp.eta_cd_eta_cnn*(newp.MUKSShabit-newp.beta_2*p.h_cd)/(1-newp.beta_2*p.h/newp.MUCSShabit)*(1-p.h/newp.MUCSShabit)/(1-p.h_cd/newp.MUKSShabit)*(1-(1-p.delta_cd)/newp.MUKSS)/newp.RCDSS;
    newp.RR=newp.eta_ch_eta_cnn*(newp.MUCSShabit-newp.beta_2*p.h_ch)/(1-newp.beta_2*p.h/newp.MUCSShabit)*(1-p.h/newp.MUCSShabit)/(1-p.h_ch/newp.MUCSShabit)*(1-(1-p.delta_ch)/newp.MUCSS)/newp.RCHSS;
    newp.Rnr=(1-(1-p.delta_)/newp.MUKSS)*newp.AA*newp.MUKSS;
    newp.ycbi_ykb=((1-p.s_AS)-newp.Rnr)/((newp.DD*(1-p.s_AS)/(1+newp.RR))+newp.Rnr);
    newp.hc_hk=newp.ycbi_ykb*(newp.RCSS*newp.MCKSS/(newp.RKSS*newp.MCCSS))^(p.alpha_/(1-p.alpha_));
    newp.HSS=0.25;
    newp.AHSS=newp.HSS;
    newp.HKSS=newp.HSS/(1+newp.hc_hk);
    newp.HCSS=newp.HSS-newp.HKSS;
    newp.HrCSS=1/3;
    newp.HrKSS=1/3;
    newp.empCSS=newp.HCSS/newp.HrCSS;
    newp.empKSS=newp.HKSS/newp.HrKSS;
    newp.ycbi=newp.HCSS*(newp.AA)^(p.alpha_/(1-p.alpha_));
    newp.ykb=newp.HKSS*(newp.AA)^(p.alpha_/(1-p.alpha_));
    newp.YCSS=newp.ycbi;
    newp.YKSS=newp.ykb;
    newp.KCSS=newp.AA*newp.ycbi*newp.MUKSS;
    newp.KKSS=newp.AA*newp.ykb*newp.MUKSS;
    newp.ECHSS=newp.RR/(1+newp.RR)*newp.ycbi*(1-p.s_AS);
    newp.ECSS=1/(1+newp.RR)*newp.ycbi*(1-p.s_AS);
    newp.ECDSS=newp.DD*newp.PKBSS*newp.ECSS;
    newp.EIKSS=(1-(1-p.delta_)/newp.MUKSS)*(newp.KCSS+newp.KKSS);
    newp.KCDSS=newp.ECDSS/(1-(1-p.delta_cd)/newp.MUKSS);
    newp.KCHSS=newp.ECHSS/(1-(1-p.delta_ch)/newp.MUCSS);
    newp.YYSS=(newp.YCSS+newp.YKSS*newp.PKBSS)/newp.PYSS;
    newp.s_k_ecd=newp.ECDSS/newp.YKSS;
    newp.s_c_ech=newp.ECHSS/newp.YCSS;
    newp.s_k_eik=newp.EIKSS/newp.YKSS;
    newp.s_yc = (newp.YCSS/newp.YYSS);
    newp.s_ecdc=newp.PKBSS*newp.ECDSS/(newp.ECSS+newp.PKBSS*newp.ECDSS+(newp.MUCSS/newp.beta_2-1+p.delta_ch)*newp.KCHSS);
    newp.INFCNASS=exp(.02/4);
    newp.INFCSS = newp.INFCNASS*((newp.MUZCSS/p.MUZKSS)^(1-p.alpha_))^(-newp.s_ecdc);
    newp.INFCORSS=newp.INFCNASS;
    newp.INFKSS=newp.INFCSS*(newp.MUZCSS/p.MUZKSS)^(1-p.alpha_);
    newp.INFWCSS=newp.INFCSS*p.MUZKSS^p.alpha_*p.MUZMSS;
    newp.INFWKSS=newp.INFWCSS;
    newp.RSS=newp.INFCSS/newp.beta_0*newp.MUCSS;
    newp.RT2SS=exp(p.tp2)*newp.RSS;
    newp.INFC10SS = newp.INFCNASS;
    newp.IMPHSSS = newp.RCHSS*newp.KCHSS;
    newp.s_k=newp.PKBSS*newp.YKSS/newp.YYSS;
    newp.INFGDPSS=newp.INFCSS^(newp.YCSS/newp.YYSS)*newp.INFKSS^(newp.YKSS*newp.PKBSS/(newp.YYSS));
    newp.LSS=newp.eta_cnn/(newp.ECSS*(1-p.h/newp.MUCSShabit))-newp.eta_cnn*newp.beta_2*p.h/(newp.ECSS*(newp.MUCSShabit-p.h));
    newp.WCSS=newp.MCCSS*(1-p.alpha_)*newp.YCSS/newp.HCSS;
    newp.WKSS=newp.MCKSS*(1-p.alpha_)*newp.YKSS/newp.HKSS;
    % xsiN_xsiH_C = ((newp.HrCSS/newp.empCSS)^(1+p.sigmah))/(1+1/p.sigmah);
    % xsiN_xsiH_K = ((newp.HrKSS/newp.empKSS)^(1+p.sigmah))/(1+1/p.sigmah);
    % gC = (1/(1+p.sigman) + 1/p.sigmah)*(xsiN_xsiH_C*(1+p.sigmah)/p.sigmah)^(-(1+p.sigman)/(1+p.sigman+p.sigmah));
    % markup_xsiN_C = (newp.HCSS^((1+p.sigmah)*(1+p.sigman)/(1+p.sigmah+p.sigman)-1))*gC/(newp.LSS*newp.WCSS);
    % gK = (1/(1+p.sigman) + 1/p.sigmah)*(xsiN_xsiH_K*(1+p.sigmah)/p.sigmah)^(-(1+p.sigman)/(1+p.sigman+p.sigmah));
    % markup_xsiN_K = (newp.HKSS^((1+p.sigmah)*(1+p.sigman)/(1+p.sigmah+p.sigman)-1))*gK/(newp.LSS*newp.WKSS);
    % markup_w = (1-newp.unempSS)^((1+p.sigmah+p.sigman)/(1+p.sigmah) - 1 - p.sigman);
    markup_w = (1-p.unempSS)^((1+p.sigmah+p.sigman)/(1+p.sigmah) - 1 - p.sigman);
    newp.theta_wc = markup_w/(markup_w -1); newp.theta_wk = newp.theta_wc;
    newp.A_HC=newp.LSS*(newp.theta_wc-1)/newp.theta_wc*newp.WCSS/(((1+p.sigman)/(1+p.sigman/(1+p.sigmah)))*newp.HCSS^(-1+(1+p.sigman)/(1+p.sigman/(1+p.sigmah))));
    newp.A_HK=newp.LSS*(newp.theta_wk-1)/newp.theta_wk*newp.WKSS/(((1+p.sigman)/(1+p.sigman/(1+p.sigmah)))*newp.HKSS^(-1+(1+p.sigman)/(1+p.sigman/(1+p.sigmah))));
    newp.xsi_NC=newp.A_HC/((1/(1+p.sigman)+1/p.sigmah)*(newp.HCSS^p.sigman/newp.HrCSS^(1+p.sigman+p.sigmah))^((1+p.sigman)/(1+p.sigman+p.sigmah)));
    newp.xsi_NK=newp.A_HK/((1/(1+p.sigman)+1/p.sigmah)*(newp.HKSS^p.sigman/newp.HrKSS^(1+p.sigman+p.sigmah))^((1+p.sigman)/(1+p.sigman+p.sigmah)));
    newp.xsi_HrC=newp.xsi_NC*(1+p.sigmah)/p.sigmah*(newp.HCSS^p.sigman/newp.HrCSS^(1+p.sigman+p.sigmah));
    newp.xsi_HrK=newp.xsi_NK*(1+p.sigmah)/p.sigmah*(newp.HKSS^p.sigman/newp.HrKSS^(1+p.sigman+p.sigmah));
    newp.UHCSS=newp.A_HC*((1+p.sigman)/(1+p.sigman/(1+p.sigmah)))*newp.HCSS^(-1+(1+p.sigman)/(1+p.sigman/(1+p.sigmah)))/newp.LSS;
    newp.UHKSS=newp.A_HK*((1+p.sigman)/(1+p.sigman/(1+p.sigmah)))*newp.HKSS^(-1+(1+p.sigman)/(1+p.sigman/(1+p.sigmah)))/newp.LSS;
    newp.HSCSS=(newp.WCSS*newp.LSS/(newp.A_HC*((1+p.sigman)/(1+p.sigman/(1+p.sigmah)))))^(1/(-1+(1+p.sigman)/(1+p.sigman/(1+p.sigmah))));
    newp.HSKSS=(newp.WKSS*newp.LSS/(newp.A_HK*((1+p.sigman)/(1+p.sigman/(1+p.sigmah)))))^(1/(-1+(1+p.sigman)/(1+p.sigman/(1+p.sigmah))));
    newp.empSCSS=((1+p.sigmah)/p.sigmah*newp.xsi_NC/newp.xsi_HrC)^(-1/(1+p.sigmah+p.sigman))*newp.HSCSS^(1/(1+p.sigman/(1+p.sigmah)));
    newp.empSKSS=((1+p.sigmah)/p.sigmah*newp.xsi_NK/newp.xsi_HrK)^(-1/(1+p.sigmah+p.sigman))*newp.HSKSS^(1/(1+p.sigman/(1+p.sigmah)));
    newp.HrSCSS=newp.HSCSS/newp.empSCSS;
    newp.HrSKSS=newp.HSKSS/newp.empSKSS;
    newp.UHSCSS=newp.A_HC*((1+p.sigman)/(1+p.sigman/(1+p.sigmah)))*newp.HSCSS^(-1+(1+p.sigman)/(1+p.sigman/(1+p.sigmah)))/newp.LSS;
    newp.UHSKSS=newp.A_HK*((1+p.sigman)/(1+p.sigman/(1+p.sigmah)))*newp.HSKSS^(-1+(1+p.sigman)/(1+p.sigman/(1+p.sigmah)))/newp.LSS;
    newp.unempSS=(newp.empSCSS+newp.empSKSS-(newp.empCSS+newp.empKSS))/(newp.empSCSS+newp.empSKSS);
    newp.QKSS=1;
    newp.QCDSS=1;
    newp.QCHSS=1;
    newp.UCSS=1;
    newp.UKSS=1;
    % XiBSS=1;
    % XiDSS=1;
    % XiHSS=1;
    newp.RL1SS=newp.RSS;
    newp.RL2SS=newp.RSS;
    newp.RL3SS=newp.RSS;
    newp.RL4SS=newp.RSS;
    newp.RL5SS=newp.RSS;
    newp.RL6SS=newp.RSS;
    newp.RL7SS=newp.RSS;
    newp.DIFFREALECSS =exp( log(newp.MUCSS));
    newp.DIFFREALEIKSS =exp( log(newp.MUKSS));
    newp.DIFFREALECDSS =exp( log(newp.MUKSS));
    newp.DIFFREALECHSS =exp( log(newp.MUCSS));
    newp.DIFFREALWSS =exp( log(newp.MUCSS) );
    newp.DIFFREALGDPSS =exp( (1-newp.s_k)*log(newp.MUCSS)+(newp.s_k)*log(newp.MUKSS));
    
    %end_steady_state;
    
    %trends;
    
    newp.DIFFREALGDPSS_obs=(1-newp.s_k)*log(newp.MUCSS)*100+(newp.s_k)*log(newp.MUKSS)*100;
    newp.DIFFREALECSS_obs=log(newp.MUCSS)*100;
    newp.DIFFREALEIKSS_obs=log(newp.MUKSS)*100;
    newp.DIFFREALECDSS_obs=log(newp.MUKSS)*100;
    newp.DIFFREALECHSS_obs=log(newp.MUCSS)*100;
    newp.DIFFREALWSS_obs=log(newp.MUCSS)*100;
    newp.INFCNASS_obs=(1-newp.s_ecdc)*log(newp.INFCSS)*100+newp.s_ecdc*log(newp.INFKSS)*100;
    newp.INFCORSS_obs=(1-newp.s_ecdc)*log(newp.INFCSS)*100+newp.s_ecdc*log(newp.INFKSS)*100;
    newp.INFKSS_obs=log(newp.INFCSS)*100-log(newp.MUKSS)*100+log(newp.MUCSS)*100;
    newp.RSS_obs=log(newp.RSS)*100;
    newp.RT2SS_obs=log(newp.RT2SS)*100;
    newp.unempSS_obs=100*log(newp.unempSS);
    
    %end_trends;
    
    ys = [
        newp.RCSS
        newp.RKSS
        newp.WCSS
        newp.WKSS
        newp.YCSS
        newp.YKSS
        newp.MCCSS
        newp.MCKSS
        newp.KCSS
        newp.KKSS
        newp.PKBSS
        newp.RSS
        newp.LSS
        newp.QKSS
        newp.HCSS
        newp.HSCSS
        newp.HKSS
        newp.HSKSS
        newp.UHCSS
        newp.UHSCSS
        newp.UHKSS
        newp.UHSKSS
        newp.empCSS
        newp.HrCSS
        newp.empKSS
        newp.HrKSS
        newp.empSCSS
        newp.HrSCSS
        newp.empSKSS
        newp.HrSKSS
        newp.unempSS
        newp.EIKSS
        newp.ECSS
        newp.INFWCSS
        newp.INFWKSS
        newp.INFCSS
        newp.INFKSS
        newp.ONE
        newp.ONE
        newp.DIFFREALGDPSS
        newp.DIFFREALECSS
        newp.DIFFREALEIKSS
        newp.DIFFREALWSS
        newp.AHSS
        newp.INFGDPSS
        newp.INFCNASS
        newp.INFCORSS
        newp.ONE
        newp.ONE
        newp.INFC10SS
        newp.ECDSS
        newp.KCDSS
        newp.RCDSS
        newp.QCDSS
        newp.KCHSS
        newp.RCHSS
        newp.ECHSS
        newp.QCHSS
        newp.KCDSS
        newp.KCHSS
        newp.USS
        newp.USS
        newp.DIFFREALECHSS
        newp.DIFFREALECDSS
        newp.beta_
        newp.ONE
        newp.ONE
        newp.ONE
        p.MUZKSS
        p.MUZMSS
        newp.ONE
        newp.MUCSS
        newp.MUKSS
        newp.ONE
        newp.ONE
        newp.ONE
        newp.RL1SS
        newp.RL2SS
        newp.RL3SS
        newp.RL4SS
        newp.RL5SS
        newp.RL6SS
        newp.RL7SS
        newp.RT2SS
        newp.DIFFREALGDPSS_obs
        newp.DIFFREALECSS_obs
        newp.DIFFREALEIKSS_obs
        newp.DIFFREALECDSS_obs
        newp.DIFFREALECHSS_obs
        newp.DIFFREALWSS_obs
        newp.ONE
        newp.INFCNASS_obs
        newp.INFCORSS_obs
        newp.INFKSS_obs
        newp.RSS_obs
        newp.RT2SS_obs
        newp.unempSS_obs
        ];
    
    
    y(id)=ys;
    
end

end