function [loglik,Incr,retcode,Filters]=switching_unscented_kalman_filter(...
syst,data_y,U,z,options)
% switching_unscented_kalman_filter - filter for nonlinear regime-swiching models
%
% ::
%
%
%   [loglik,Incr,retcode,Filters]=switching_unscented_kalman_filter(...
%    syst,data_y,U,z,options)
%
% Args:
%
%    - **syst** [struct]: structure containing:
%
%          - **PAI00** [vector]: initial probability distributions of regimes
%
%          - **a** [cell]: initial conditions in each regime
%
%          - **Qfunc** [function handle]: transition matrix generator
%
%          - **ff** [function handle]: ft=ff(rt,xt,et), where rt is the
%          regime, xt is the vector of state variables and et the vector of
%          shocks
%
%          - **P** [cell]: initial covariance matrix of the states in each
%          regime
%
%          - **H** [cell]: Measurement error covariance matrices in each regime
%
%          - **SIGeta** [cell]: Covariance matrix of structural shocks.
%
%    - **data_y** [matrix]: ny x T matrix of data
%
%    - **U** [[]|matrix]: ndx x T matrix of exogenous data
%
%    - **z** [function handle|logical|vector]: linear connection of the
%    observables to the state.
%
%    - **include_in_likelihood** [logical]: selector of increments to include
%    in the likelihood calculation (temporarily disabled)
%
%    - **options** [struct]: structure with various options
%
% Returns:
%    :
%
%    - **loglik** [scalar]: log likelihood
%
%    - **Incr** [vector]: increments of elements going into the likelihood
%
%    - **retcode** [{0}|integer]: flag for problems.
%
%    - **Filters** [struct]: Filtered, updated and smoothed variables
%
% Note:
%
% Example:
%
%    See also:

% data
%-----
% obs_id=syst.obs_id;

% first iterate to consider in the computation of the likelihood.
%-------------------------------------------------------------------
first=syst.start;

% state matrices
%---------------
ff=syst.ff;

% T=syst.Tx;

R=syst.Te;

H=syst.H;

Qfunc=syst.Qfunc;
% sep_compl=syst.sep_compl;
% cond_shocks_id=syst.anticipated_shocks;
% ss=syst.steady_state;
% xlocs=syst.state_vars_location;

% initial conditions
%-------------------
a=syst.a;

P=syst.P;

PAItt=syst.PAI00;
% RR=syst.Te_Te_prime;
% current size of the system after expansion
m=syst.m; % size(T{1},1);
% size of the system prior to expansion
% m_orig=syst.m_orig;

% free up memory
%---------------
clear syst

h=numel(a);

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

nshocks=exo_nbr*horizon;

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

% definitions and options
%------------------------
twopi=2*pi;

% smoothing not ready
%---------------------
store_filters=options.kf_filtering_level;

if store_filters
    
    nsteps=options.kf_nsteps;
    
    if store_filters>2
        
        % place holder for smoothing gains
        %----------------------------------
        SG=repmat({nan(m,m,smpl)},1,h);
        
    end
    
else
    % do not do multi-step forecasting during estimation
    nsteps=1;
    
end

kalman_tol=options.kf_tol;

% initialization of matrices
%-----------------------------
loglik=[];

Incr=nan(smpl,1);

[Filters]=utils.filtering.initialize_storage(a,P,PAI,Q,p0,exo_nbr,horizon,nsteps,smpl,store_filters);

twopi_p_dF=nan(1,h);
% the following elements change size depending on whether observations are
% missing or not and so it is better to have them in cells rather than
% matrices
iPvv=cell(1,h);

v=cell(1,h);

% This also changes size but we need to assess whether we reach the steady state fast or not
K=zeros(m,p0,h); % <---K=cell(1,h);

% no problem
%-----------
retcode=0;

% lump states and shocks together
%---------------------------------
mm=m+nshocks;

k=3-mm;

w=1/(2*(mm+k))*ones(1,2*mm+1);

w(1)=k/(mm+k);

sqrt_m_plus_k=sqrt(mm+k);

for t=1:smpl% <-- t=0; while t<smpl,t=t+1;
    % data and indices for observed variables at time t
    %--------------------------------------------------
    [p,occur,obsOccur]=z(t);
    
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
        
        Pav=P{st}(:,obsOccur); % PZt=<-- P{st}*Z';
        
        Pvv=Pav(obsOccur,:); % <--- F=Z*PZt+H{st}(occur,occur);
        
        if h_last>0
            
            Pvv=Pvv+H{st}(occur,occur,min(t,h_last));
            
        end
        
        detF=det(Pvv);
        
        failed=detF<=0;
        
        if ~failed
            
            iPvv{st}=Pvv\eye(p);
            
            failed=any(isnan(iPvv{st}(:)));
            
        end
        
        if failed
            
            retcode=305;
            
            return
            
        end
        
        % Kalman gain (for update)
        %-------------------------
        K(:,occur,st)=Pav*iPvv{st}; % K=PZt/F{st};
        
        % state covariance update (Ptt=P-P*Z*inv(F)*Z'*P)
        %------------------------------------------------
        P{st}=P{st}-K(:,occur,st)*Pav.';%<---P{st}=P{st}-K(:,occur,st)*P{st}(obsOccur,:);
        
        twopi_p_dF(st)=twopi^p*detF;
        
        % state update (att=a+K*v)
        %-------------------------
        a{st}=a{st}+K(:,occur,st)*v{st};
        
        % likelihood contribution
        %-------------------------
        log_f01(st)=-0.5*(...
            log(twopi_p_dF(st))+...
            v{st}'*iPvv{st}*v{st}...
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
    
    % Likelihood computation
    %-----------------------
    Incr(t)=log(likt);
    
    % endogenous probabilities (conditional on time t information)
    %-------------------------------------------------------------
    att=a;
    
    Ptt=P;
    
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
        
        P{splus}=zeros(m);
        
        for st=1:h
            
            if h==1
                
                pai_st_splus=1;
                
            else
                
                pai_st_splus=Q(st,splus)*PAItt(st)/PAI(splus);
                
            end
            
            a{splus}=a{splus}+pai_st_splus*att{st};
            
            P{splus}=P{splus}+pai_st_splus*Ptt{st};
            
        end
        
        [xtt,retcode]=utils.filtering.sigma_points(a{splus},P{splus},'ukf',sqrt_m_plus_k);
		
        if retcode
            
			return
            
        end
        
        a_plus=0;
        
        x_plus=[xtt,zeros(m,2*nshocks)];
        
        % run sigma points for states 
        %-----------------------------
        for jj=1:2*m+1
            
            [x_plus(:,jj),~,retcode]=ff(splus,xtt(:,jj),shocks,Ut);
            
            if retcode
                
                return
                
            end
            
            a_plus=a_plus+w(jj)*x_plus(:,jj);
            
        end
        
        % run sigma points for shocks: we keep the state itself at a{splus}
        %------------------------------------------------------------------
        offset=2*m+1;
        
        for jj=1:nshocks
            
            this_shocks=shocks;
            
            this_shocks(jj)=1*sqrt_m_plus_k;
            
            [x_plus(:,offset+1),~,retcode]=ff(splus,a{splus},this_shocks,Ut);
            
            if retcode
                
                return
                
            end
            
            a_plus=a_plus+w(offset+1)*x_plus(:,offset+1);
            
            offset=offset+1;
            
            this_shocks(jj)=-1*sqrt_m_plus_k;
            
            [x_plus(:,offset+1),~,retcode]=ff(splus,a{splus},this_shocks,Ut);
            
            if retcode
                
                return
                
            end
            
            a_plus=a_plus+w(offset+1)*x_plus(:,offset+1);
            
            offset=offset+1;
            
        end
        
        P_plus=0;
        
        C_t_tp1=0;
        
        for jj=1:size(x_plus,2)
            
            x_a=x_plus(:,jj)-a_plus;
            
            P_plus=P_plus+w(jj)*(x_a*x_a.');
            
            if store_filters>2 && jj<=2*m+1
                
                C_t_tp1=C_t_tp1+w(jj)*(xtt(:,jj)-a{splus})*x_a.';
                % the remaining increments are zeros since the first
                % parenthesis in these cases is a{splus}-a{splus}
                
            end
            
        end
        
        a{splus}=a_plus;
        
        % the symmetrization/projection step could have been avoided but we
        % want to give a chance to the covariance matrix of forecast errors
        % to be well-behaved.
        %------------------------------------------------------------------
        P_plus=utils.cov.project(P_plus);
        
        P{splus}=P_plus;
        
        if store_filters>2
            
            % Compute smoothing gain and store for later use...
            SG{splus}(:,:,t)=C_t_tp1*pinv(P_plus);
            
        end
        
    end
    
    if store_filters>0
        
        store_predictions()
        
    end
    
end

% included only if in range
loglik=sum(Incr(first:end));

if store_filters>2 % store smoothed
    
    for t=smpl:-1:1
        
        Q=Filters.Q(:,:,t);
        
        [~,occur,obsOccur]=z(t);
        
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
            
            % smoothed state
            %------------------
            Filters.atT{s0}(:,1,t)=Filters.att{s0}(:,1,t);
            % P_{t|T}=P_{t|t}+SG*(P_{t+1|T}-P_{t+1|t})*SG.'
            
            if t<smpl
                
                Filters.atT{s0}(:,1,t)=Filters.atT{s0}(:,1,t)+...
                    SG{s0}(:,:,t)*(Filters.atT{s0}(:,1,t+1)-Filters.a{s0}(:,1,t+1));
            
            end
            
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


    function store_updates()
        
        Filters.PAItt(:,t)=PAItt;
        
        for st_=1:h
            
            Filters.att{st_}(:,1,t)=a{st_};
            
            Filters.Ptt{st_}(:,:,t)=P{st_};
%             if store_filters>2
%                 K_store{st_}(:,occur,t)=K(:,occur,st_);
%                 iF_store{st_}(occur,occur,t)=iPvv{st_};
%                 v_store{st_}(occur,t)=v{st_};
%             end
        
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