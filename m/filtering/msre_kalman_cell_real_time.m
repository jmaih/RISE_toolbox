function [loglik,Incr,retcode,Filters]=msre_kalman_cell_real_time(syst,data_info,state_trend,init,options)
% H1 line
%
% Syntax
% -------
% ::
%
% Inputs
% -------
%
% Outputs
% --------
%
% More About
% ------------
%
% Examples
% ---------
%
% See also:


% this filter assumes a state space of the form
% X_t=c_t{st}+T{st}*X_{t-1}+R{st}*eps_t
% y_t=d_t{st}+Z*X_t+eta_t
% where c_t{st} and d_t{st} are, possibly time-varying, deterministic terms
% the covariance matrices can be time-varying

% data
%-----
data_structure=data_info.data_structure;
include_in_likelihood=data_info.include_in_likelihood;
no_more_missing=data_info.no_more_missing;
obs_id=data_info.varobs_id;
data=data_info.y;
dataz=data_info.z;
last_good_conditional_observation=data_info.last_good_conditional_observation;
% N.B: data_info also contains x, the observations on the exogenous
% variables (trend, etc). But those observations will come through
% data_trend and so are not used directly in the filtering function.

% state matrices
%---------------
T=syst.T;
R=syst.R;
H=syst.H;
Qfunc=syst.Qfunc;

% initial conditions
%-------------------
a=init.a;
P=init.P;
PAItt=init.PAI00;
% RR=init.RR;
h=numel(T);
a_nsteps=cell(1,h);

% free up memory
%---------------
m=size(T{1},1);
npages=data_info.npages;
[~,exo_nbr,horizon]=size(R{1});
kmax=horizon-1;
horizon=min(npages-1,horizon); % the first page is hard information
ncp=horizon;
hypothesis=options.forecast_conditional_hypothesis;
kf_nan_future_obs_means_missing=options.kf_nan_future_obs_means_missing;
state_vars_location=1:m;
        
% holder of conditional data
%----------------------------
y0=struct('y',nan(m,1,horizon+1));

restr_y_id_in_y=imag(data_info.restr_y_id);
restr_y_id_in_state=real(data_info.restr_y_id);
restr_z_id_in_state=real(data_info.restr_z_id);
% if ~isempty(restr_z_id_in_state)
%     error('conditioning on shocks is not allowed in estimation')
% end
% n_y_rest=numel(restr_y_id_in_y);

DPHI=cell(1,h);
DT=cell(1,h);
ss=cell(1,h);
Im=eye(m);
for st=1:h
    [DPHI{st},DT{st}]=utils.forecast.conditional.build_shock_restrictions(...
        T{st},R{st},...
        restr_y_id_in_state,restr_z_id_in_state,...
        ncp,... data availability
        horizon,... horizon +or-1
        hypothesis);
    % redress: remove the third dimension
    %-------------------------------------
    R{st}=R{st}(:,:);
    if st==1
        nc=size(DT{st},2);
        nshocs=ncp*numel(restr_z_id_in_state);
    end
    DT{st}=[DT{st};zeros(nshocs,nc)];
    ss{st}=(Im-T{st})\state_trend{st}(:,1);
end

Q=Qfunc(a{1});
PAI=transpose(Q)*PAItt;

% matrices' sizes
%----------------
[p0,smpl]=size(data);
smpl=min(smpl,find(include_in_likelihood,1,'last'));

h_last=0;
if ~isempty(H{1})
    h_last=size(H{1},3);
end

% few transformations
%--------------------
Tt=T;
any_T=false(1,m);
for st=1:h
    Tt{st}=transpose(T{st}); % permute(T,[2,1,3]); % transpose
    any_T=any_T|any(abs(T{st})>1e-9,1);
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
    'simul_update_shocks_handle',[],...
    'simul_do_update_shocks',false,...
    'forecast_conditional_hypothesis',options.forecast_conditional_hypothesis,...
    'nsteps',nsteps,...
    'burn',0,...
    'k_future',kmax,...
    'y',[]);

kalman_tol=options.kf_tol;
expanded_flag=store_filters>2;
if expanded_flag
    shock_span=size(DPHI{1},2);
    mm=m+shock_span;
    % remember all shocks have standard deviation 1
    %----------------------------------------------
    Pstandard=eye(mm);
    for st=1:h
        a{st}=[a{st};zeros(shock_span,1)];
        Pstandard(1:m,1:m)=P{st};
        P{st}=Pstandard;
    end
    clear Pstandard
else
    mm=m;
end
clear data_info syst init

% initialization of matrices
%-----------------------------
loglik=[];
Incr=nan(smpl,1);

if expanded_flag
    K_store=[];
    iF_store=[];
    v_store=[];
    Tt_store=repmat({nan(mm,mm,smpl)},1,h);
    Rt_store=repmat({nan(mm,shock_span,smpl)},1,h);
end
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
K=zeros(mm,p0,h); % <---K=cell(1,h);

% no problem
%-----------
retcode=0;
is_steady=false;

% disp(['do not forget to test whether it is possible to reach the steady ',...
%     'state with markov switching'])

for t=1:smpl% <-- t=0; while t<smpl,t=t+1;
    % data and indices for observed variables at time t
    %--------------------------------------------------
    occur=data_structure(:,t);
    p=sum(occur); % number of observables to be used in likelihood computation
    y=data(occur,t);
    obsOccur=obs_id(occur); %<-- Z=M.Z(occur,:);
    
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
    
    [Q,retcode]=Qfunc(att{1});
    if retcode
        return
    end
    
    % Probabilities predictions
    %--------------------------
    if h>1
        PAI=Q'*PAItt;
    end
    
    % advance information: skip the first page with hard information
    %----------------------------------------------------------------
    MUt=data(restr_y_id_in_y,t,1+(1:horizon));
    shocks_MUt=dataz(:,t,1+(1:horizon));
    
    lgcobs=min(last_good_conditional_observation(t),horizon);
    % keep the contemporaneous by default
    lgcobs=max(lgcobs,1);

    % state and state covariance prediction
    %--------------------------------------
    for splus=1:h
        [a_nsteps{splus},P{splus},Tt{splus},Rt{splus}]=predict_while_collapsing(expanded_flag);
        % trim in case many forecasts where produced in earlier rounds
        %--------------------------------------------------------------
        a{st}(1:m)=a_nsteps{st}(:,1);
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
loglik=sum(Incr(include_in_likelihood));

if expanded_flag % store smoothed
    r=zeros(mm,h);
    ZZ=eye(mm);
    ZZ=ZZ(obs_id,:);
    for t=smpl:-1:1
        Q=Filters.Q(:,:,t);
        occur=data_structure(:,t);
        obsOccur=obs_id(occur); %<-- Z=M.Z(occur,:);
        Z=ZZ(occur,:);
        y=data(occur,t);
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

% let's play squash
%-------------------
if store_filters>0
    for st=1:h
        %         Filters.filt_shocks{st}=Filters.a{st}(m+1:end,:,:);
        Filters.a{st}=Filters.a{st}(1:m,:,:); % second dimension is the real-time forecasting steps
        Filters.P{st}=Filters.P{st}(1:m,1:m,:);
        if store_filters>1
            %             Filters.update_shocks{st}=Filters.att{st}(m+1:end,:,:);
            Filters.att{st}=Filters.att{st}(1:m,:,:); % second dimension is the real-time forecasting steps
            Filters.Ptt{st}=Filters.Ptt{st}(1:m,1:m,:);
            if store_filters>2
                % In real-time, the smoothed shocks of interest are to be
                % found after the endogenous variables in the state vector
                % and not in the usual smoothed disturbances. This is
                % because in the presence of conditional information,
                % smoothed shocks are no longer IID whereas smoothed
                % disturbances are. In order to be able to reconstruct the
                % smoothed series in the absence of the time varying state
                % matrices, we need the smoothed shocks, not the smoothed
                % disturbances. And so we discard the smoothed disturbances
                Filters.eta{st}=Filters.atT{st}(m+1:end,:,:);
                Filters.atT{st}=Filters.atT{st}(1:m,:,:); % second dimension is the real-time forecasting steps
                Filters.PtT{st}=Filters.PtT{st}(1:m,1:m,:);
            end
        end
    end
end

    function [a,P,Tt,Rt]=predict_while_collapsing(ExpandedFlag)
        a=0;
        P=0;
        Tt=0;
        Rt=0;
        for snow=1:h
            pai_snow_Over_slead=Q(snow,splus)*PAItt(snow)/PAI(splus);
            [a00_,P00_,Tt_snow,Rt_snow]=prediction_step(...
                T{splus},R{splus},att{snow},Ptt{snow},...
                OMGt,DPHI{splus},DT{splus},ExpandedFlag);
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
            if expanded_flag
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
            Filters.a{splus_}(1:m,:,t+1)=a_nsteps{splus_};
            Filters.P{splus_}(:,:,t+1)=P{splus_};
        end
    end
    function Filters=initialize_storage()
        Filters=struct();
        if store_filters>0 % store filters
            Filters.a=repmat({zeros(mm,nsteps,smpl+1)},1,h);
            Filters.P=repmat({zeros(mm,mm,smpl+1)},1,h);
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
                Filters.att=repmat({zeros(mm,1,smpl)},1,h);
                Filters.Ptt=repmat({zeros(mm,mm,smpl)},1,h);
                Filters.PAItt=zeros(h,smpl);
                if expanded_flag % store smoothed
                    K_store=repmat({zeros(mm,p0,smpl)},1,h);
                    iF_store=repmat({zeros(p0,p0,smpl)},1,h);
                    v_store=repmat({zeros(p0,smpl)},1,h);
                    Filters.atT=repmat({zeros(mm,1,smpl)},1,h);
                    Filters.PtT=repmat({zeros(mm,mm,smpl)},1,h);
                    Filters.eta=repmat({zeros(shock_span,1,smpl)},1,h); % smoothed shocks
                    Filters.epsilon=repmat({zeros(p0,1,smpl)},1,h); % smoothed measurement errors
                    Filters.PAItT=zeros(h,smpl);
                end
            end
        end
    end
    function [a,P,Tt,Rt]=prediction_step(T,R,att,Ptt,OMGt,DPHI,DT,ExpandedFlag)
        narginchk(4,8);
        reduced=nargin==4;
        % zero the useless entries otherwise the variables will appear as missing
        if ~kf_nan_future_obs_means_missing
            R(1:m,lgcobs*exo_nbr+1:end)=0;
        end
        Tt=T;
        Rt=R;
        bt=0;
        MUt_splus=MUt(:,:);
        MUt_splus=MUt_splus-ss{splus}(restr_y_id_in_state)*ones(1,horizon);
        MUt_splus_old=[MUt_splus(:)
            shocks_MUt(:)];
        if ~reduced
            reduced=reduced||(~ExpandedFlag && all(isnan(MUt_splus_old(:))));
        end
        
        if reduced
            % no valid future information available, nothing to match
            Record=[];
        else
            [Tt,Rt,bt,~,Record]=utils.forecast.conditional.state_matrices(T,R,MUt_splus_old,OMGt,DPHI,DT,Record,ExpandedFlag);
        end
        
        RR=Rt*Rt';
        P=Tt*Ptt*transpose(Tt)+RR;
        % Make sure we remain symmetric
        P=utils.cov.symmetrize(P);
        
        if nsteps==1
            a=bt+Tt*att;
        else
            shocks(restr_z_id_in_state,:)=shocks_MUt(:,:);
            fkst_options_.shocks=shocks;
            
            y0.y(:,1,1)=att(state_vars_location);
            y0.y(restr_y_id_in_state,1,2:end)=MUt(:,:);
            [a,~,retcode]=utils.forecast.multi_step(y0,ss(splus),...
                {[T,zeros(m,1),R]},state_vars_location,fkst_options_);
        end
    end
end