function [loglik,Incr,retcode,Filters]=constrained_regime_switching_kalman_filter_cell(...
syst,data_y,U,z,options)%syst,data_info,state_trend,init,options


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
ff=syst.ff;

T=syst.Tx;

R=syst.Te;

H=syst.H;

Qfunc=syst.Qfunc;

sep_compl=syst.sep_compl;

% cond_shocks_id=syst.anticipated_shocks;
% ss=syst.steady_state;
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
m_orig=syst.m_orig;

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

% matrices' sizes
%----------------
[p0,smpl]=size(data_y);

[~,exo_nbr,horizon]=size(R{1});

shocks=zeros(exo_nbr,horizon);

% deterministic variables only if rows present
%----------------------------------------------
tmax_u=size(U,2)*(size(U,1)>0);

if tmax_u
    
    % update the state (a_{1|0} with the deterministic variables
    %------------------------------------------------------------
    Ut=U(:,0+1);
    
    for splus=1:h
        
        [a{splus}]=ff(splus,a{splus},shocks,Ut);
        
    end
    
end

h_last=0;

if ~isempty(H{1})
    
    h_last=size(H{1},3);
    
end

rqr_last=size(RR{1},3);

if rqr_last>1
    
    error('time-varying impact matrices not supported in this filter')
    
end

% definitions and options
%------------------------
twopi=2*pi;

store_filters=options.kf_filtering_level;

if store_filters
    
    nsteps=options.kf_nsteps;
    
    if store_filters>2
        
        R_store=struct();
        
    end
     
else
    
    % do not do multi-step forecasting during estimation
    nsteps=1;
    
end

kalman_tol=options.kf_tol;

% time-varying R matrix
%-----------------------
Rt=cell(1,h);

% few transformations
%--------------------
Tt=T;

for st=1:h
    
    Tt{st}=transpose(T{st}); % permute(T,[2,1,3]); % transpose
    
end

% initialization of matrices
%-----------------------------
loglik=[];

Incr=nan(smpl,1);

[Filters,K_store,iF_store,v_store]=utils.filtering.initialize_storage(a,P,PAI,Q,p0,exo_nbr,horizon,nsteps,smpl,store_filters);

oldK=inf;

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
    [p,occur,obsOccur,no_more_missing]=z(t);
    
    y=data_y(occur,t);
    
    log_f01 = nan(h,1);
    
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
            P{st}=P{st}-K(:,occur,st)*PZt.';%<---P{st}=P{st}-K(:,occur,st)*P{st}(obsOccur,:);
            
            twopi_p_dF(st)=twopi^p*detF;
            
        end
        % state update (att=a+K*v)
        %-------------------------
        a{st}=a{st}+K(:,occur,st)*v{st};
        
        log_f01(st)=-0.5*(...
            log(twopi_p_dF(st))+...
            v{st}'*iF{st}*v{st}...
            );
        
    end
    
    [Incr(t),PAI01_tt,retcode]=switch_like_exp_facility(PAI,log_f01,kalman_tol);
    
    if retcode
        
        return
        
    end
    
    PAItt=sum(PAI01_tt,2);
    
    if store_filters>1
        
        store_updates();
        
    end
    
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
    
    % state and state covariance prediction (this is next period and so,
    % the exogenous variables should be picked for period t+1
    %--------------------------------------------------------------------
    if t+1<=tmax_u
        Ut=U(:,t+1);
    else
        Ut=[];
    end
    for splus=1:h
        
        a{splus}=zeros(m,1);
        
        if ~is_steady
            
            P{splus}=zeros(m);
            
        end
        
        for st=1:h
            
            if h==1
                
                pai_st_splus=1;
                
            else
                
                pai_st_splus=Q(st,splus)*PAItt(st)/PAI(splus);
                
            end
            
            a{splus}=a{splus}+pai_st_splus*att{st};
            
            if ~is_steady
                
                P{splus}=P{splus}+pai_st_splus*Ptt{st};
                
            end
            
        end
        
        [a{splus},is_active_shock,retcode]=ff(splus,a{splus},shocks,Ut);
        
        if retcode
            
            return
            
        end
        
        Rt{splus}=R{splus};
        
        if ~is_steady
            % we will never be steady if we have restrictions!!!
            %---------------------------------------------------
            if isempty(sep_compl)
                
                RR_splus=RR{splus};
                
            else
                % always force the first shock to be active
                is_active_shock(1)=true;
                
                % Now kill the inactive locations
                Rt{splus}(1:m_orig,:,~is_active_shock)=0;
                
                RR_splus=Rt{splus}(:,:)*Rt{splus}(:,:).';
                
            end
            
            P{splus}=T{splus}*P{splus}(xlocs,xlocs)*Tt{splus}+RR_splus;
            
            P{splus}=utils.cov.symmetrize(P{splus});
            
        end
        
    end
    
    if store_filters>0
        
        store_predictions()
        
        if store_filters>2
            
            R_store(t).R=Rt;
            
        end
        
    end
    
    if ~is_steady && isempty(sep_compl) % && h==1
        % don't check steadiness if the filter is constrained.
        %------------------------------------------------------
        
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
                K_store{s0}(:,occur,t),Filters.P{s0}(:,:,t),resquare(T{s0}),...
                R_store(t).R{s0}(:,:),Z,iF_store{s0}(occur,occur,t),...
                v_store{s0}(occur,t));
            
            % smoothed measurement errors
            %--------------------------
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

    function T=resquare(T)
        
        tmp=zeros(m);
        
        tmp(:,xlocs)=T;
        
        T=tmp;
        
    end

    function store_updates()
        
        Filters.PAItt(:,t)=PAItt;
        
        for st_=1:h
            
            Filters.att{st_}(:,1,t)=a{st_};
            
            Filters.Ptt{st_}(:,:,t)=P{st_};
            
            if store_filters>2
                
                K_store{st_}(:,occur,t)=K(:,occur,st_);
                
                iF_store{st_}(occur,occur,t)=iF{st_};
                
                v_store{st_}(occur,t)=v{st_};
                
            end
            
        end
        
    end

    function store_predictions()
        
        Filters.PAI(:,t+1)=PAI;
        
        Filters.Q(:,:,t+1)=Q;
        
        for splus_=1:h
            
            Filters.a{splus_}(:,1,t+1)=a{splus_};
            
            Filters.P{splus_}(:,:,t+1)=P{splus_};
            
            for istep_=2:nsteps
                % this assumes that we stay in the same state and we know
                % we will stay. The more general case where we can jump to
                % another state is left to the forecasting routine.
                if t+istep_-1<=tmax_u
                    
                    Utplus=U(:,t+istep_-1);
                    
                else
                    
                    Utplus=[];
                    
                end
                
                [Filters.a{splus_}(:,istep_,t+1),~,rcode]=ff(splus_,...
                    Filters.a{splus_}(:,istep_-1,t+1),shocks,Utplus);
                
                if rcode
                    % do not exit completely...
                    break
                    
                end
                
            end
            
        end
        
    end

end