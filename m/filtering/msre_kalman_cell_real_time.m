function [loglik,Incr,retcode,Filters]=msre_kalman_cell_real_time(...
syst,y_data,U,z,e_data,options)

% H1 line
%
% ::
%
%
% Args:
%
% Returns:
%    :
%
% Note:
%
% Example:
%
%    See also:


% this filter assumes a state space of the form
% X_t=c_t{st}+T{st}*X_{t-1}+R{st}*eps_t
% y_t=d_t{st}+Z*X_t+eta_t
% where c_t{st} and d_t{st} are, possibly time-varying, deterministic terms
% the covariance matrices can be time-varying

% data
%-----
obs_id=syst.obs_id;

% first iterate to consider in the computation of the likelihood.
%-------------------------------------------------------------------
first=syst.start;

% state matrices
%---------------
% ff=syst.ff;

T=syst.Tx;

R=syst.Te;

H=syst.H;

Qfunc=syst.Qfunc;
% sep_compl=syst.sep_compl;
xlocs=syst.state_vars_location;

% initial conditions
%-------------------
a=syst.a;

P=syst.P;

PAItt=syst.PAI00;

RR=syst.Te_Te_prime;
% current size of the system after expansion
m=syst.m; % size(T{1},1);
% size of the system prior to expansion
% m_orig=syst.m_orig;
ss=syst.steady_state;

% free up memory
%---------------
clear syst

h=numel(T);

att_all=a{1}*PAItt(1);

for rt=2:h
    
    att_all=att_all+a{rt}*PAItt(rt);
    
end

Q=Qfunc(att_all);

PAI=transpose(Q)*PAItt;

% data
%------
data_y=y_data{1};

dataz=e_data{1};

restr_y_id_in_y=imag(y_data{2});

restr_y_id_in_state=real(y_data{2});

restr_z_id_in_state=real(e_data{2});

% matrices' sizes
%----------------
[p0,smpl,npages]=size(data_y);

[~,exo_nbr,horizon]=size(R{1});

kmax=horizon-1;

horizon=min(npages-1,horizon); % the first page is hard information

ncp=horizon;

hypothesis=options.forecast_conditional_hypothesis;

if ~strcmpi(hypothesis,'nas')

    error('real-time filtering no longer works with hypotheses other than NAS')

end

kf_nan_future_obs_means_missing=options.kf_nan_future_obs_means_missing;

[nrows,tmax_u]=size(U);

if min(nrows,tmax_u)>0

    error('real time filtering with exogenous variables not ready')

end

% this system works with square matrices for the time being
%-----------------------------------------------------------
tmp=zeros(m);

for ireg=1:h

    tmp(:,xlocs)=T{ireg};
   
    T{ireg}=tmp;

end

clear tmp

% initial conditions
%-------------------
a_nsteps=cell(1,h);

% free up memory
%---------------
state_vars_location=1:m;
        
% holder of conditional data
%----------------------------
y0=struct('y',nan(m,1,horizon+1));

active_shocks=setdiff(1:exo_nbr,restr_z_id_in_state);

DPHI=cell(1,h);

DT=cell(1,h);

% Conditon on all the existing shocks and nan the relevant ones where
% necessary

for st=1:h

    [DPHI{st},DT{st}]=utils.filtering.build_shock_restrictions(...
        T{st},R{st},...
        restr_y_id_in_state,1:exo_nbr,... restr_z_id_in_state
        ncp);% ... data availability
        
    if st==1
    
        nc=size(DT{st},2);
        
        nshocks=ncp*exo_nbr;%<--nshocs=ncp*numel(restr_z_id_in_state);
    
    end
    % redress: remove the third dimension
    %-------------------------------------
    R{st}=R{st}(:,:);
    
    DT{st}=[DT{st};zeros(nshocks,nc)];

end

h_last=0;

if ~isempty(H{1})

    h_last=size(H{1},3);

end

% Intialize time-varying matrices
%---------------------------------
Tt=T;

Rt=R;

OMGt=[]; % uncertainty of the conditions

Record=[];

% definitions and options
%------------------------
twopi=2*pi;

store_filters=options.kf_filtering_level;

if store_filters

    nsteps=options.kf_nsteps;

else
    % do not do multi-step forecasting during estimation
    nsteps=1;

end

shocks=nan(exo_nbr,horizon);

fkst_options_=struct('PAI',1,...
    'Qfunc',@(x)1,...
    'states',1,...
    'shocks',shocks,...
    'forecast_conditional_hypothesis',options.forecast_conditional_hypothesis,...
    'nsteps',nsteps,...
    'burn',0,...
    'k_future',kmax,...
    'y',[]);

kalman_tol=options.kf_tol;

if store_filters

    shock_span=size(DPHI{1},2);
    
    K_store=[];
    
    iF_store=[];
    
    v_store=[];
    
    Tt_store=repmat({nan(m,m,smpl)},1,h);
    
    Rt_store=repmat({nan(m,shock_span,smpl)},1,h);

end

% initialization of matrices
%-----------------------------
loglik=[];

Incr=nan(smpl,1);

Filters=initialize_storage();

oldK=inf;

PAI01y=nan(h,1);

twopi_p_dF=nan(1,h);

% the following elements change size depending on whether observations are
% missing or not and so it is better to have them in cells rather than
% matrices
iF=cell(1,h);

v=cell(1,h);
% This also changes size but we need to assess whether we reach the steady state fast or not

K=zeros(m,p0,h); % <---K=cell(1,h);

% no problem
%-----------
retcode=0;

is_steady=false;

% disp(['do not forget to test whether it is possible to reach the steady ',...
%     'state with markov switching'])

for t=1:smpl% <-- t=0; while t<smpl,t=t+1;
    % data and indices for observed variables at time t
    %--------------------------------------------------
    [p,occur,obsOccur,no_more_missing,lgcobs_t]=z(t);

    y=data_y(occur,t,1);
    
    likt=0;
    
    for st=1:h
        % forecast of observables: already include information about the
        % trend and or the steady state from initialization
        %------------------------------------------------------------------
        yf=a{st}(obsOccur); %<-- yf=Z*a{st};
        
        % forecast errors and variance
        %-----------------------------
        v{st}=y-yf;
        
        if ~is_steady
    
            PZt=P{st}(:,obsOccur); % PZt=<-- P{st}*Z';
            
            Fst=PZt(obsOccur,:); % <--- F=Z*PZt+H{st}(occur,occur);
            
            if h_last>0
            
                Fst=Fst+H{st}(occur,occur,min(t,h_last));
            
            end
            
            detF=det(Fst);
            
            failed=detF<=0;
            
            if ~failed
            
                iF{st}=Fst\eye(p);
                
                failed=any(isnan(iF{st}(:)));
            
            end
            
            if failed
            
                retcode=305;
                
                return
            
            end
            
            % Kalman gain (for update)
            %-------------------------
            K(:,occur,st)=PZt*iF{st}; % K=PZt/F{st};
            
            % state covariance update (Ptt=P-P*Z*inv(F)*Z'*P)
            %------------------------------------------------
            P{st}=P{st}-K(:,occur,st)*PZt.';
            
            twopi_p_dF(st)=twopi^p*detF;
        
        end
        % state update (att=a+K*v)
        %-------------------------
        a{st}=a{st}+K(:,occur,st)*v{st};
        
        f01=(twopi_p_dF(st)*exp(v{st}'*iF{st}*v{st}))^(-0.5);
        
        PAI01y(st)=PAI(st)*f01;
        
        likt=likt+PAI01y(st);
    
    end
    
    % Probability updates
    %--------------------
    PAI01_tt=PAI01y/likt;
    
    if likt<kalman_tol && (any(isnan(PAI01_tt))||any(isinf(PAI01_tt)))
    
        retcode=306;
        
        return
    
    end
    
    PAItt=sum(PAI01_tt,2);
    
    % Likelihood computation
    %-----------------------
    Incr(t)=log(likt);
    
    % endogenous probabilities (conditional on time t information)
    %-------------------------------------------------------------
    att=a;
    
    if ~is_steady
    
        Ptt=P;
    
    end
    
    if h>1
        
        att_all=att{1}*PAItt(1);
        
        for rt=2:h
            
            att_all=att_all+att{rt}*PAItt(rt);
            
        end
        
        [Q,retcode]=Qfunc(att_all);
        
        if retcode
        
            return
        
        end
        
        % Probabilities predictions
        %--------------------------
        PAI=Q'*PAItt;
        
    end
    
    % advance information: skip the first page with hard information
    %----------------------------------------------------------------
    MUt=data_y(restr_y_id_in_y,t,1+(1:horizon));
    
    MUt=MUt(:,:);

    shocks_MUt=dataz(:,t,1+(1:horizon));
    
    last_mu=find(any(~isnan(MUt),1),1,'last');
    
    if isempty(last_mu)
        
        last_mu=0;
        
    end
    
    lgcobs=min(lgcobs_t,horizon);
    
    lgcobs=min(last_mu,lgcobs);
    
    % keep the contemporaneous by default
    lgcobs=max(lgcobs,1);

    % state and state covariance prediction
    %--------------------------------------
    for splus=1:h
        
        [a_nsteps{splus},P{splus},Tt{splus},Rt{splus}]=predict_while_collapsing();
        
        % trim in case many forecasts where produced in earlier rounds
        %--------------------------------------------------------------
        a{st}=a_nsteps{st}(:,1);
        
    end
    
    if store_filters>0
        
        if store_filters>1
            
            store_updates();
            
        end
        
        store_predictions()
        
    end
    
    if ~is_steady % && h==1
        
        [is_steady,oldK]=utils.filtering.check_steady_state_kalman(...
            is_steady,K,oldK,options,t,no_more_missing);
        
    end
    
end

% included only if in range
loglik=sum(Incr(first:end));

if store_filters>2 % store smoothed
    
    r=zeros(m,h);
    
    ZZ=eye(m);
    
    ZZ=ZZ(obs_id,:);
    
    for t=smpl:-1:1
        
        Q=Filters.Q(:,:,t);
        
        [~,occur,obsOccur]=z(t);
        
        Z=ZZ(occur,:);
        
        y=data_y(occur,t);
        
        for s0=1:h
            
            for s1=1:h
                % joint probability of s0 (today) and s1 (tomorrow)
                if t==smpl
                    
                    pai_0_1=Q(s0,s1)*Filters.PAItt(s0,t);
                    
                else
                    
                    pai_0_1=Q(s0,s1)*Filters.PAItt(s0,t)*...
                        Filters.PAItT(s1,t+1)/Filters.PAI(s1,t+1);
                    
                end
                
                % smoothed probabilities
                %-----------------------
                Filters.PAItT(s0,t)=Filters.PAItT(s0,t)+pai_0_1;
                
            end
            
            % smoothed state and shocks
            %--------------------------
            [Filters.atT{s0}(:,1,t),Filters.eta{s0}(:,1,t),r(:,s0)]=...
                utils.filtering.smoothing_step(Filters.a{s0}(:,1,t),r(:,s0),...
                K_store{s0}(:,occur,t),Filters.P{s0}(:,:,t),...
                Tt_store{s0}(:,:,t),Rt_store{s0}(:,:,t),...% T{s0},R{s0},
                Z,iF_store{s0}(occur,occur,t),v_store{s0}(occur,t));
            
            % smoothed measurement errors
            %------------------------------
            Filters.epsilon{s0}(occur,1,t)=y-Filters.atT{s0}(obsOccur,1,t);
            
        end
        % correction for the smoothed probabilities [the approximation involved does not always work
        % especially when dealing with endogenous switching.
        SumProbs=sum(Filters.PAItT(:,t));
        
        if abs(SumProbs-1)>1e-8
            
            Filters.PAItT(:,t)=Filters.PAItT(:,t)/SumProbs;
            
        end
        
    end
    
end

    function [a,P,Tt,Rt]=predict_while_collapsing()
        
        a=0;
        
        P=0;
        
        Tt=0;
        
        Rt=0;
        
        for snow=1:h
            
            pai_snow_Over_slead=Q(snow,splus)*PAItt(snow)/PAI(splus);
            [a00_,P00_,Tt_snow,Rt_snow]=prediction_step(...
                T{splus},R{splus},att{snow},Ptt{snow},...
                OMGt,DPHI{splus},DT{splus});
            
            a=a+pai_snow_Over_slead*a00_;
            
            P=P+pai_snow_Over_slead*P00_;
            
            Tt=Tt+pai_snow_Over_slead*Tt_snow;
            
            Rt=Rt+pai_snow_Over_slead*Rt_snow;
            
        end
        
    end

    function store_updates()
        
        Filters.PAItt(:,t)=PAItt;
        
        for st_=1:h
            
            Filters.att{st_}(:,1,t)=att{st_};
            
            Filters.Ptt{st_}(:,:,t)=Ptt{st_};
            
            if store_filters
                
                K_store{st_}(:,occur,t)=K(:,occur,st_);
                
                iF_store{st_}(occur,occur,t)=iF{st_};
                
                v_store{st_}(occur,t)=v{st_};
                
                % state matrices
                Tt_store{st_}(:,:,t)=Tt{st_};
                
                Rt_store{st_}(:,:,t)=Rt{st};
                
            end
            
        end
        
    end

    function store_predictions()
        
        Filters.PAI(:,t+1)=PAI;
        
        Filters.Q(:,:,t+1)=Q;
        
        for splus_=1:h
            
            Filters.a{splus_}(:,:,t+1)=a_nsteps{splus_};
            
            Filters.P{splus_}(:,:,t+1)=P{splus_};
            
        end
        
    end

    function Filters=initialize_storage()
        
        Filters=struct();
        
        if store_filters>0 % store filters
            
            Filters.eta_tlag=cell(1,h);
            
            Filters.a=repmat({zeros(m,nsteps,smpl+1)},1,h);
            
            Filters.P=repmat({zeros(m,m,smpl+1)},1,h);
            
            for state=1:h
                
                Filters.a{state}(:,1,1)=a{state};
                
                Filters.P{state}(:,:,1)=P{state};
                
            end
            
            Filters.PAI=zeros(h,smpl+1);
            
            Filters.PAI(:,1)=PAI;
            
            for istep=2:nsteps
                
                % in steady state, we remain at the steady state
                %------------------------------------------------
                for state=1:h
                    
                    Filters.a{state}(:,istep,1)=Filters.a{state}(:,istep-1,1);
                    
                end
                
            end
            
            Filters.Q=zeros(h,h,smpl+1);
            
            Filters.Q(:,:,1)=Q;
            
            if store_filters>1 % store updates
                
                Filters.eta_tt=cell(1,h);
                
                Filters.att=repmat({zeros(m,1,smpl)},1,h);
                
                Filters.Ptt=repmat({zeros(m,m,smpl)},1,h);
                
                Filters.PAItt=zeros(h,smpl);
                
                if store_filters>2 % store smoothed
                    
                    K_store=repmat({zeros(m,p0,smpl)},1,h);
                    
                    iF_store=repmat({zeros(p0,p0,smpl)},1,h);
                    
                    v_store=repmat({zeros(p0,smpl)},1,h);
                    
                    Filters.atT=repmat({zeros(m,1,smpl)},1,h);
                    
                    Filters.PtT=repmat({zeros(m,m,smpl)},1,h);
                    
                    Filters.eta=repmat({zeros(shock_span,1,smpl)},1,h); % smoothed shocks
                    
                    Filters.epsilon=repmat({zeros(p0,1,smpl)},1,h); % smoothed measurement errors
                    
                    Filters.PAItT=zeros(h,smpl);
                    
                end
                
            end
            
        end
        
    end

    function [a,P,Tt,Rt]=prediction_step(T,R,att,Ptt,OMGt,DPHI,DT)
        
        narginchk(7,7);
        
        % zero the useless entries otherwise the variables will appear as missing
        shocks_=nan(exo_nbr,horizon);
        
        shocks_(restr_z_id_in_state,:)=shocks_MUt(:,:);
        
        if ~kf_nan_future_obs_means_missing
            
            shocks_(active_shocks,lgcobs+1:end)=0;
            
            R(1:m,lgcobs*exo_nbr+1:end)=0;
            
        end
        
        MUt_splus=MUt-ss{splus}(restr_y_id_in_state)*ones(1,horizon);
        
        MUt_splus_old=[MUt_splus(:)
            shocks_(:)];
        
        ExpandedFlag=false; 
        
        [Tt,Rt,bt,~,Record]=utils.forecast.rscond.state_matrices(T,R,MUt_splus_old,OMGt,DPHI,DT,Record,ExpandedFlag);
        
        RR=Rt*Rt';
        
        P=Tt*Ptt*transpose(Tt)+RR;
        
        % Make sure we remain symmetric
        P=utils.cov.symmetrize(P);
        
        % do the first step
            a=ss{splus}+bt+Tt*(att-ss{splus});
            
        % add the subsequent steps artificially in order to match the sizes
        if nsteps>1
            
            fkst_options_.shocks=shocks_;
            
            y0.y(:,1,1)=att(state_vars_location);
            
            y0.y(restr_y_id_in_state,1,2:end)=MUt(:,:);
            
            [a_multi_steps,~,retcode]=utils.forecast.multi_step(y0,ss(splus),...
                {[T,zeros(m,1),R]},state_vars_location,fkst_options_);
            
            a=a(:,ones(1,nsteps));
            
            a(1:m,:)=a_multi_steps;
            
            % future shocks not collected for multi-step forecasting
            a(m+1:end,2:end)=nan; 
            
        end
        
    end

end