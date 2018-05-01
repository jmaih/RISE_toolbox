function [loglik,Incr,retcode,Filters]=msre_kalman(syst,data_info,state_trend,init,options)
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
data_structure=data_info.data_structure;

include_in_likelihood=data_info.include_in_likelihood;

no_more_missing=data_info.no_more_missing;

obs_id=data_info.varobs_id;

data=data_info.y;
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

RR=init.RR;

% free up memory
%---------------
clear data_info syst init

[m,m,h]=size(T);

att_all=a(:,1)*PAItt(1);

for rt=2:h
    
    att_all=att_all+a(:,rt)*PAItt(rt);
    
end

Q=Qfunc(att_all);

PAI=transpose(Q)*PAItt;

% matrices' sizes
%----------------
[p0,smpl]=size(data);

smpl=min(smpl,find(include_in_likelihood,1,'last'));

nshocks=size(R,2);

c_last=0;if ~isempty(state_trend),c_last=size(state_trend,2);end

rqr_last=size(RR,4);

h_last=0;

if ~isempty(H)
    
    h_last=size(H,4);
    
end

% few transformations
%--------------------
Tt=permute(T,[2,1,3]);

any_T=false(1,m);

for st=1:h
    
    any_T=any_T|any(abs(T(:,:,st))>1e-9,1);
    
end

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

kalman_tol=options.kf_tol;

% initialization of matrices
%-----------------------------
loglik=[];

Incr=nan(smpl,1);

if store_filters>2
    
    K_store=[];
    
    iF_store=[];
    
    v_store=[];
    
end

Filters=initialize_storage();

oldK=inf;

twopi_p_dF=nan(1,h);

K=zeros(m,p0,h); % <---K=cell(1,h);

% no problem
%-----------
retcode=0;

is_steady=false;

for t=1:smpl% <-- t=0; while t<smpl,t=t+1;
    % data and indices for observed variables at time t
    %--------------------------------------------------
    occur=data_structure(:,t);
    
    p=sum(occur); % number of observables to be used in likelihood computation
    
    y=data(occur,t);
    
    obsOccur=obs_id(occur); %<-- Z=M.Z(occur,:);
    
    if store_filters>2
        
        v=nan(p,h);
        
        iF=nan(p,p,h);
        
    end
    
    log_f01 = nan(h,1);
    
    for st=1:h
        % forecast of observables: already include information about the
        % trend and or the steady state from initialization
        %------------------------------------------------------------------
        yf=a(obsOccur,st); %<-- yf=Z*a{st};
        
        % forecast errors and variance
        %-----------------------------
        v_st=y-yf;
        
        if ~is_steady
            
            PZt=P(:,obsOccur,st); % PZt=<-- P{st}*Z';
            
            Fst=PZt(obsOccur,:); % <--- F=Z*PZt+H{st}(occur,occur);
            
            if h_last>0
                
                Fst=Fst+H(occur,occur,st,min(t,h_last));
                
            end
            
            detF=det(Fst);
            
            failed=detF<=0;
            
            if ~failed
                
                iF_st=Fst\eye(p);
                
                failed=any(isnan(iF_st(:)));
                
            end
            
            if failed
                
                retcode=305;
                
                return
                
            end
            
            % Kalman gain (for update)
            %-------------------------
            K(:,occur,st)=PZt*iF_st; % K=PZt/F{st};
            
            % state covariance update (Ptt=P-P*Z*inv(F)*Z'*P)
            %------------------------------------------------
            P(:,:,st)=P(:,:,st)-K(:,occur,st)*PZt.';
            
            twopi_p_dF(st)=twopi^p*detF;
            
        end
        % state update (att=a+K*v)
        %-------------------------
        a(:,st)=a(:,st)+K(:,occur,st)*v_st;
        
        log_f01(st)=-0.5*(...
            log(twopi_p_dF(st))+...
            v{st}'*iF{st}*v{st}...
            );
        
        if store_filters>2
            
            iF(:,:,st)=iF_st;
            
            v(:,st)=v_st;
            
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
        
        att_all=att(:,1)*PAItt(1);
        
        for rt=2:h
            
            att_all=att_all+att(:,rt)*PAItt(rt);
            
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
    a(:)=0;%<---a(:,:)=0;
    
    if ~is_steady
        
        P(:,:,:)=0;
        
    end
    
    for splus=1:h
        
        for st=1:h
            
            if h==1
                
                pai_st_splus=1;
                
            else
                
                pai_st_splus=Q(st,splus)*PAItt(st)/PAI(splus);
                
            end
            
            a(:,splus)=a(:,splus)+pai_st_splus*att(:,st);
            
            if ~is_steady
                
                P(:,:,splus)=P(:,:,splus)+pai_st_splus*Ptt(:,:,st);
            
            end
            
        end
        
        a(:,splus)=T(:,any_T,splus)*a(any_T,splus); % a{splus}=T{splus}*a{splus};
        
        if ~is_steady
            
            P(:,:,splus)=T(:,:,splus)*P(:,:,splus)*Tt(:,:,splus)+RR(:,:,splus,min(t,rqr_last));
            %             P(:,:,splus)=T(:,any_T,splus)*P(any_T,any_T,splus)*Tt(any_T,:,splus)+RR(:,:,splus,min(t,rqr_last));
            P(:,:,splus)=utils.cov.symmetrize(P(:,:,splus));
        
        end
        
        if c_last>0 % <-- ~isempty(state_trend)
            
            a(:,splus)=a(:,splus)+state_trend(:,min(t+1,c_last),splus);
            
        end
        
    end
    
    if store_filters>0
        
        store_predictions()
        
    end
    
    if ~is_steady % && h==1
        
        [is_steady,oldK]=utils.filtering.check_steady_state_kalman(...
            is_steady,K,oldK,options,t,no_more_missing);
        
    end
    
end

% included only if in range
loglik=sum(Incr(include_in_likelihood));

if store_filters>2 % store smoothed
    
    r=zeros(m,h);
    
    ZZ=eye(m);
    
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
                K_store{s0}(:,occur,t),Filters.P{s0}(:,:,t),T(:,:,s0),R(:,:,s0),Z,...
                iF_store{s0}(occur,occur,t),v_store{s0}(occur,t));
            
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
            
            Filters.att{st_}(:,1,t)=a(:,st_);
            
            Filters.Ptt{st_}(:,:,t)=P(:,:,st_);
            
            if store_filters>2
                
                K_store{st_}(:,occur,t)=K(:,occur,st_);
                
                iF_store{st_}(occur,occur,t)=iF(:,:,st_);
                
                v_store{st_}(occur,t)=v(:,st_);
                
            end
            
        end
        
    end

    function store_predictions()
        
        Filters.PAI(:,t+1)=PAI;
        
        Filters.Q(:,:,t+1)=Q;
        
        for splus_=1:h
            
            Filters.a{splus_}(:,1,t+1)=a(:,splus_);
            
            Filters.P{splus_}(:,:,t+1)=P(:,:,splus_);
            
            for istep_=2:nsteps
                % this assumes that we stay in the same state and we know
                % we will stay. The more general case where we can jump to
                % another state is left to the forecasting routine.
                Filters.a{splus_}(:,istep_,t+1)=T(:,:,splus_)*Filters.a{splus_}(:,istep_-1,t+1);
                
                if c_last>0 % <-- ~isempty(state_trend)
                    
                    Filters.a{splus_}(:,istep_,t+1)=Filters.a{splus_}(:,istep_,t+1)+state_trend(:,min(t+istep_,c_last),splus_);
                
                end
                
            end
            
        end
        
    end

    function Filters=initialize_storage()
        
        Filters=struct();
        
        if store_filters>0 % store filters
            
            Filters.a=repmat({zeros(m,nsteps,smpl+1)},1,h);
            
            Filters.P=repmat({zeros(m,m,smpl+1)},1,h);
            
            for state=1:h
                
                Filters.a{state}(:,1,1)=a(:,state);
                
                Filters.P{state}(:,:,1)=P(:,:,state);
                
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
                
                Filters.att=repmat({zeros(m,1,smpl)},1,h);
                
                Filters.Ptt=repmat({zeros(m,m,smpl)},1,h);
                
                Filters.PAItt=zeros(h,smpl);
                
                if store_filters>2 % store smoothed
                    
                    K_store=repmat({zeros(m,p0,smpl)},1,h);
                    
                    iF_store=repmat({zeros(p0,p0,smpl)},1,h);
                    
                    v_store=repmat({zeros(p0,smpl)},1,h);
                    
                    Filters.atT=repmat({zeros(m,1,smpl)},1,h);
                    
                    Filters.PtT=repmat({zeros(m,m,smpl)},1,h);
                    
                    Filters.eta=repmat({zeros(nshocks,1,smpl)},1,h); % smoothed shocks
                    
                    Filters.epsilon=repmat({zeros(p0,1,smpl)},1,h); % smoothed measurement errors
                    
                    Filters.PAItT=zeros(h,smpl);
                    
                end
                
            end
            
        end
        
    end

end