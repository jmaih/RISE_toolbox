function [loglik,Incr,retcode,Filters]=crs_linear_filter_update_cell_estim_K(...
syst,data_y,U,z,options,impose_conditions)

% crs_linear_filter_update_cell_estim_K -- filter with update of K. gain
%
% ::
%
%
%   [loglik,Incr,retcode,Filters]=crs_linear_filter_update_cell_estim_K(...
%    syst,y,U,z,options)
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
%    - **y** [matrix]: ny x T x npages matrix of data
%
%    - **U** [[]|matrix]: ndx x T matrix of exogenous data
%
%    - **z** [function handle|logical|vector]: linear connection of the
%    observables to the state.
%
%    - **include_in_likelihood** [logical]: selector of increments to include
%    in the likelihood calculation
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
%    - The filter checks for violations of constraints and uses the
%    information in the database to reset the offending variables to their
%    values in the database. When this occurs, the Kalman gain is recomputed
%    so that the traditional updating equation holds.
%
%    - This strategy is adopted so as to permit the use of the efficient
%    smoothing algorithm of Durbin and Koopman instead of the classical
%    smoothing algorithm which requires multiple inversions of a potentially
%    singular covariance matrix.
%
% Example:
%
%    See also:


% this filter assumes a state space of the form
% X_t=c_t{st}+T{st}*X_{t-1}+R{st}*eps_t
% y_t=d_t{st}+Z*X_t+eta_t
% where c_t{st} and d_t{st} are, possibly time-varying, deterministic terms
% the covariance matrices can be time-varying

cutoff=-sqrt(eps);

if nargin<6
    
    impose_conditions=false;
    
end
% data
%-----
obs_id=syst.obs_id;

% first iterate to consider in the computation of the likelihood.
%-------------------------------------------------------------------
first=syst.start;

% regime-wise anticipation
%--------------------------
k=syst.k(:).';

has_fire=k>0;

% state matrices
%---------------
ff=syst.ff;

T=syst.Tx;

Tbig=syst.T;

R=syst.Te;

H=syst.H;

Qfunc=syst.Qfunc;

sep_compl=syst.sep_compl;

cond_shocks_id=syst.anticipated_shocks;

ss=syst.steady_state;

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
[p0,smpl,npages]=size(data_y);

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

myshocks=cell(1,h);

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

% this will be useful for multi-step forecasting
%------------------------------------------------
expected_shocks=cell(1,h);

% no problem
%-----------
retcode=0;

is_steady=false;

% update the options for conditional forecasting
%------------------------------------------------
options.PAI=1;

options.states=ones(horizon,1);

options.Qfunc=Qfunc;

options.y=[];

options.burn=0;

options.k_future=horizon-1;

% forecast all the way...
options.nsteps=horizon;

% originally, the update is the same as the filter since a=T*a . If this is
% not the case, it does not matter much since this is just the beginning of
% the filter.
%--------------------------------------------------------------------------
att=a;

for t=1:smpl% <-- t=0; while t<smpl,t=t+1;
    % data and indices for observed variables at time t
    %--------------------------------------------------
    [p,occur,obsOccur,no_more_missing]=z(t);
    
    y=data_y(occur,t,1);
    
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
        if impose_conditions
            
            [a{st},K(:,occur,st),retcode,myshocks{st}]=state_update_without_test(a{st},K(:,occur,st),v{st},st);
        
        else
            
            [a{st},K(:,occur,st),retcode,myshocks{st}]=state_update_with_test(a{st},K(:,occur,st),v{st},st);
            % <--- a{st}=a{st}+K(:,occur,st)*v{st};
        end
        
        if retcode
            
            return
            
        end
        
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
    
    % state and state covariance prediction
    %--------------------------------------
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
        
        expected_shocks{splus}=shocks;
        
        for st=1:h
            
            if h==1
                
                pai_st_splus=1;
                
            else
                
                pai_st_splus=Q(st,splus)*PAItt(st);
                
                if pai_st_splus>0
                    
                    pai_st_splus=pai_st_splus/PAI(splus);
                    
                end
                
            end
            
            a{splus}=a{splus}+pai_st_splus*att{st};
            % forecast with the expected shocks we had from the
            % updating step
            %---------------------------------------------------
            expected_shocks{splus}=expected_shocks{splus}+...
                pai_st_splus*myshocks{st}(:,1:horizon);
            
            if st==h
                % remove the first period since the shocks were
                % expected already from the period before...
                expected_shocks{splus}=[expected_shocks{splus}(:,2:horizon),shocks(:,1)];
            
            end
            
            if ~is_steady
                
                P{splus}=P{splus}+pai_st_splus*Ptt{st};
                
            end
            
        end
        
        [a{splus},is_active_shock,retcode]=ff(splus,a{splus},expected_shocks{splus},Ut);
        
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

    function y0=simul_initial_conditions(a_filt,start_iter)
        
        if nargin<2
            
            start_iter=horizon;
            
        end
        
        options.shocks=shocks;
        
        % all shocks can be used in the update (first period)
        %----------------------------------------------------
        options.shocks(:,1)=nan;
        
        % beyond the first period, only the shocks with long reach can
        % be used
        options.shocks(cond_shocks_id,1:start_iter)=nan;
        
        % compute a conditional forecast
        ycond=reshape(data_y(:,t,1:start_iter),p,[]);
        
        ycond=struct('data',ycond(:,:,ones(3,1)),'pos',obs_id);
        
        econd=[options.shocks,zeros(exo_nbr,start_iter-1)];
        
        econd=struct('data',econd(:,:,ones(3,1)),'pos',1:exo_nbr);
        
        rcond=struct('data',ones(start_iter,1),'pos',nan);
        
        y0=struct('y',a_filt,'ycond',ycond,'econd',econd,'rcond',rcond);
        
    end

    function [a_update,K,retcode,myshocks]=state_update_without_test(a_filt,K,v,st)
        
        retcode=0;
        
        lgcobs=last_good_future_information(t);
        
        if lgcobs==0
            % compute the simple update: expected shocks are zero
            %-----------------------------------------------------
            atmp=a_filt+K*v;
            
            myshocks=shocks;
            
            is_redo_K=false;
            
        else
            
            is_redo_K=true;
            % compute the conditional update: note that we are now using
            % the last update as the initial condition
            %-------------------------------------------------------------
            y0=simul_initial_conditions(att{st});
            
            options_=options;
            
            options_.nsteps=lgcobs+1;
            
            options_.shocks(:,options_.nsteps+1:end)=0;
            
            [fsteps,~,retcode,~,myshocks]=utils.forecast.multi_step(y0,...
                ss(st),Tbig(st),xlocs,options_);
            
            atmp=fsteps(:,1);
            
        end
        
        a_update=atmp;
        
        if is_redo_K
            
            K=recompute_kalman_gain(a_update,a_filt,v);
            
        end
        
        function lgcobs=last_good_future_information(t)
            
            expected_data=squeeze(data_y(:,t,2:min(horizon,npages)));
            
            good_obs=any(~isnan(expected_data),1);
            
            lgcobs=length(good_obs);
            
            while ~isempty(good_obs)
                
                if good_obs(end)
                    
                    break
                    
                else
                    
                    good_obs(end)=[];
                    
                    lgcobs=lgcobs-1;
                    
                end
                
            end
            
        end
        
    end

    function [a_update,retcode,myshocks_,violations]=state_update_with_test(a_filt,Kv,st)
        
        retcode=0;
        % compute the update
        %--------------------
        atmp=a_filt+Kv;
        
        % compute one-step forecast from the update and check whether we
        % can expect violations
        %----------------------------------------------------------------
        atmp_ss=atmp-ss{st};
        
        a_expect=ss{st}+T{st}*atmp_ss(xlocs);
        
        violations=~isempty(sep_compl) && any(sep_compl(a_expect)<cutoff);
        
        % initialize the shocks
        %------------------------
        myshocks_=shocks;
        
        a_update=atmp;
        
        if has_fire(st) && violations
            
            do_anticipation()
            
        end
        
        function do_anticipation()
            % if we can expect violations, inform the state that constraints
            % will be binding
            %----------------------------------------------------------------
            if options.debug
                
                old_update=atmp(:,ones(1,horizon));
                
                old_update(:,2:end)=nan;
                
            end
            
            start_iter=1;
            
            while start_iter<horizon && violations
                
                start_iter=start_iter+1;
                % compute the conditional update: note that we are now
                % using the last update as the initial condition
                %----------------------------------------------------------  
                
                y0=simul_initial_conditions(att{st},start_iter);
                
                myoptions=options;
                
                myoptions.nsteps=start_iter;
                
                myoptions.states=options.states(1:start_iter);
                
                [fsteps,~,retcode,~,myshocks_]=utils.forecast.multi_step(y0,ss(st),Tbig(st),xlocs,myoptions);
                
                if retcode
                    
                    a_update=[];
                    
                    return
                    
                end
                
                for iter=1:myoptions.nsteps
                    
                    violations=any(sep_compl(fsteps(:,iter))<cutoff);
                    
                    if violations
                        
                        break
                        
                    end
                    
                end
                
                if start_iter==horizon
                    % if we have come so far, then there is no violation in the
                    % last step. But there could be some in the future step
                    f__=fsteps(:,end)-ss{st};
                    
                    a_expect=ss{st}+T{st}*f__(xlocs);
                    
                    violations=any(sep_compl(a_expect)<cutoff);
                
                end
                
                if ~violations
                    
                    atmp=fsteps(:,1);
                    
                end
                
                if options.debug
                    
                    old_update(:,start_iter)=fsteps(:,1);
                    
                end
                
            end
            % make sure we have the correct size since there may be more shocks
            % resulting from forecasting multi-periods
            %------------------------------------------------------------------
            myshocks_=myshocks_(:,1:horizon);
            
            if options.debug && start_iter>1
                
                keyboard
                
            end
            
            if violations
                
                retcode=701;
                
            end
            
            a_update=atmp;
            
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
            
            % remove the first period since the shocks were
            % expected already from the period before...
            shocks_splus=[expected_shocks{splus}(:,2:horizon),shocks(:,1)];
            
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
                    Filters.a{splus_}(:,istep_-1,t+1),shocks_splus,Utplus);
                
                % update the shocks
                shocks_splus=[shocks_splus(:,2:end),shocks(:,1)];
                
                if rcode
                    % do not exit completely...
                    break
                    
                end
                
            end
            
        end
        
    end

end

function K=recompute_kalman_gain(a_update,a_filt,v)

K=(a_update-a_filt)*pinv(v);

end