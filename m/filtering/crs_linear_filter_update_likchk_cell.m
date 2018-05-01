function [loglik,Incr,retcode,Filters]=crs_linear_filter_update_likchk_cell(...
syst,data_y,U,z,options,impose_conditions)


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

is_recompute_K=false;

quick_collapse=true;

check_shock_viol=false; % check the violation of the shocks sign

check_is_violation_in_period_2=true;

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
    
    error('deterministic variables currently disabled')
    
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

[Filters]=utils.filtering.initialize_storage(a,P,PAI,Q,p0,exo_nbr,horizon,nsteps,smpl,store_filters);

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

a_expect=cell(1,h);

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

is_violation=false(smpl,h);
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
        
        ast_old=a{st};
        
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
            
            [a{st},a_expect{st},retcode,myshocks{st},K(:,occur,st)]=state_update_without_test(a{st},K(:,occur,st),v{st},st);
        
        else
            
            [a{st},a_expect{st},retcode,myshocks{st},K(:,occur,st)]=state_update_with_test(a{st},K(:,occur,st),v{st},st);
            % <--- a{st}=a{st}+K(:,occur,st)*v{st};
        end
        
        if retcode
            
            return
            
        end
        
        if check_is_violation_in_period_2
            
            aviol=a_expect{st};
            
        else
            
            aviol=ast_old;
            
        end
        
        is_violation(t,st)=~isempty(sep_compl) && any(sep_compl(aviol)<cutoff);
        
        if is_violation(t,st)
            
            log_f01(st)=-inf;
            
        else
            
            log_f01(st)=-0.5*(...
                log(twopi_p_dF(st))+...
                v{st}'*iF{st}*v{st}...
                );
            
        end
        
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
        
        shk=0;
        
        for st=1:h
            
            if h==1
                
                pai_st_splus=1;
                
            else
                
                pai_st_splus=Q(st,splus)*PAItt(st);
                
                if pai_st_splus>0
                    
                    pai_st_splus=pai_st_splus/PAI(splus);
                    
                end
                
            end
            
            if quick_collapse
                
                a{splus}=a{splus}+pai_st_splus*a_expect{st}(:,1);
                
                shk=shk+pai_st_splus*myshocks{st};
                
            else
                
                a{splus}=a{splus}+pai_st_splus*att{st};
                
            end
            
            if ~is_steady
                
                P{splus}=P{splus}+pai_st_splus*Ptt{st};
                
            end
            
        end
        
        if quick_collapse
            
            retcode=0;
            
            expected_shocks{splus}=[shk(:,2:end),zeros(exo_nbr,1)];
            
            is_active_shock=any(abs(expected_shocks{splus})>sqrt(eps),1);
            
        else
            
            [a{splus},expected_shocks{splus},is_active_shock,retcode]=compute_onestep(a{splus},splus);
            
            if retcode
                
                return
                
            end
            
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
            % smoothed state and shocks
            %--------------------------
            if t==smpl
                
                Filters.atT{s0}(:,1,t)=Filters.att{s0}(:,1,t);
                
            else
                
                [Filters.atT{s0}(:,1,t)]=utils.filtering.smoothing_step_classical(...
                    resquare(T{s0}),Filters.att{s0}(:,1,t),...
                    Filters.Ptt{s0}(:,:,t),Filters.P{s0}(:,:,t+1),...
                    Filters.atT{s0}(:,1,t+1),Filters.a{s0}(:,1,t+1));
                
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

    function [a_update,a_expect,retcode,myshocks,K]=state_update_without_test(a_filt,K,v,st)
        
        retcode=0;
        
        lgcobs=last_good_future_information(t);
        
        if lgcobs==0
            % compute the simple update: expected shocks are zero
            %-----------------------------------------------------
            atmp=a_filt+K*v;
            
            myshocks=shocks;
            
            atmp_ss=atmp-ss{st};
            
            a_expect=ss{st}+T{st}*atmp_ss(xlocs);
            
        else
            % compute the conditional update: note that we are now using
            % the last update as the initial condition
            %-------------------------------------------------------------
            y0=crs_filter_simul_init(att{st},shocks,cond_shocks_id,...
                data_y,t,p,obs_id,true);
            
            options_=options;
            
            options_.nsteps=lgcobs+1;
            
            options_.shocks(:,options_.nsteps+1:end)=0;
            
            [fsteps,~,retcode,~,myshocks]=utils.forecast.multi_step(y0,...
                ss(st),Tbig(st),xlocs,options_);
            
            atmp=fsteps(:,1);
            
            a_expect=fsteps(:,2);
            
            if is_recompute_K
                
                K=recompute_kalman_gain(atmp,a_filt,v);
                
            end
            
        end
        
        a_update=atmp;
        
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

    function [a_update,fkst,retcode,myshocks_,K,violations]=state_update_with_test(a_filt,K,v,st)
        % compute the update
        %--------------------
        atmp=a_filt+K*v;
        % compute one-step forecast from the update and check whether we
        % can expect violations
        %----------------------------------------------------------------
        atmp_ss=atmp-ss{st};
        
        a_expect_=ss{st}+T{st}*atmp_ss(xlocs);
        
        violations=~isempty(sep_compl) && any(sep_compl(a_expect_)<cutoff);
        % initialize the shocks
        %------------------------
        myshocks_=shocks;
        
        retcode=0;
        
        a_update=atmp;
        
        if has_fire(st) && violations
            
            nsteps__=2;
            
            match_first_page=true;
            % compute the conditional update: note that we are now
            % using the last update as the initial condition
            %----------------------------------------------------------
            [a_update,a_expect_,myshocks_,retcode]=do_anticipation(att{st},nsteps__,...
                match_first_page);
            
            if is_recompute_K
                
                K=recompute_kalman_gain(a_update,a_filt,v);
                
            end
            
        end
        
        if retcode
            
            fkst=[];
            
        else
            
            fkst=roll_forward(a_update,[myshocks_(:,2:end),zeros(exo_nbr,1)],nsteps,st);
        
        end
        
    end

    function fkst=roll_forward(a0,shocks,nsteps,st)
        
        fkst=a0(:,ones(nsteps,1));
        
        for istep=1:nsteps
            
            a0_ss=a0-ss{st};
            
            a0=ss{st}+T{st}*a0_ss(xlocs);
            
            a0=a0+R{st}(:,:)*shocks(:);
            
            fkst(:,istep)=a0;
            
            shocks=[shocks(:,2:end),zeros(exo_nbr,1)];
            
        end
        
    end

    function [a1,myshocks_,is_active_shock,retcode]=compute_onestep(a0,st,...
            check_second_too)
        
        if nargin<3
            
            check_second_too=false;
            
        end
        % compute one-step forecast from the update and check whether we
        % can expect violations
        %----------------------------------------------------------------
        a0_ss=a0-ss{st};
        
        a1=ss{st}+T{st}*a0_ss(xlocs);
        
        violations=~isempty(sep_compl) && any(sep_compl(a1)<cutoff);
        
        if check_second_too
            
            a1_ss=a1-ss{st};
            
            a2=ss{st}+T{st}*a1_ss(xlocs);
            
            violations=violations||...
                (~isempty(sep_compl) && any(sep_compl(a2)<cutoff));
            
        end
        
        retcode=0;
        
        myshocks_=shocks;
        
        % initialize default of is_active_shock
        is_active_shock=any(abs(shocks)>sqrt(eps),1);
        
        if has_fire(st) && violations
            
            nsteps__=1+check_second_too;
            
            match_first_page=false;
            
            [a1,~,myshocks_,retcode]=do_anticipation(a0,nsteps__,match_first_page);
            
            is_active_shock=any(abs(myshocks_)>sqrt(eps),1);
            
        end
        
    end

    function [a1,a_expect,myshocks_,retcode]=do_anticipation(a0,nsteps,match_first_page)
        
        check_first_shock=match_first_page;
        retcode=0;
        
        a1=[];
        
        a_expect=[];
        
        myshocks_=[];
        
        if options.debug
            
            record_forecasts=nan*a0;
            
            record_forecasts=record_forecasts(:,ones(1,horizon),ones(1,horizon));
       
        end
        % if we can expect violations, inform the state that constraints
        % will be binding
        %----------------------------------------------------------------
        is_viol=true;
        
        layers=0;
        
        nsteps=nsteps-1;
        
        if check_first_shock
            
            start_shocks=1;
            
        else
            
            start_shocks=2;
            
        end
        
        while nsteps<horizon && is_viol
            
            nsteps=nsteps+1;
            
            layers=layers+1;
            
            y0=crs_filter_simul_init(a0,shocks,cond_shocks_id,data_y,t,p,...
                obs_id,match_first_page,nsteps);
            
            myoptions=options;
            
            myoptions.nsteps=nsteps;
            
            myoptions.states=options.states(1:nsteps);
            
            [fsteps,~,retcode,~,myshocks_]=utils.forecast.multi_step(y0,ss(st),Tbig(st),xlocs,myoptions);
            
            if retcode
                
                return
                
            end
            
            is_viol=false;
            
            if check_shock_viol
                
                is_viol=any(vec(myshocks_(cond_shocks_id,start_shocks:end))<cutoff);
            
            end
            
            if ~is_viol
                
                for iter=1:myoptions.nsteps
                    
                    is_viol=any(sep_compl(fsteps(:,iter))<cutoff);
                    
                    if is_viol
                        
                        break
                        
                    end
                    
                end
                
                if ~is_viol
                    % if we have come so far, then there is no violation in the
                    % last step. But there could be some in the future step
                    f__=fsteps(:,end)-ss{st};
                    
                    a_extra=ss{st}+T{st}*f__(xlocs);
                    
                    is_viol=any(sep_compl(a_extra)<cutoff);
                    
                end
                
            end
            
            if ~is_viol
                
                a1=fsteps(:,1);
                
                if nsteps>1
                    
                    a_expect=fsteps(:,2);
                    
                else
                    
                    a_expect=a_extra;
                    
                end
                
            end
            
            if options.debug
                
                record_forecasts(:,:,layers)=[fsteps,nan(nrows,horizon-nsteps)];
                
                if nsteps==horizon
                    
                    record_forecasts=record_forecasts(:,:,1:layers);
                    
                end
                
            end
            
        end
        % make sure we have the correct size since there may be more shocks
        % resulting from forecasting multi-periods
        %------------------------------------------------------------------
        myshocks_=myshocks_(:,1:horizon);
        
        if options.debug && nsteps>1
            
            keyboard
            
        end
        
        if is_viol
            
            retcode=701;
            
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
            
        end
        
    end

    function store_predictions()
        
        Filters.PAI(:,t+1)=PAI;
        
        Filters.Q(:,:,t+1)=Q;
        
        for splus_=1:h
            
            if quick_collapse
                
                Filters.a{splus_}(:,1:nsteps,t+1)=a_expect{splus_};
                
                continue
                
            else
                
                Filters.a{splus_}(:,1,t+1)=a{splus_};
                
                Filters.P{splus_}(:,:,t+1)=P{splus_};
            
            end
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
