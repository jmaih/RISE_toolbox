function [loglik,Incr,retcode,Filters]=msre_kalman_cell(syst,data_info,data_trend,state_trend,init,options)

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
Q=syst.Q0;

% initial conditions
%-------------------
a=init.a;
P=init.P;
PAItt=init.PAI00;
RR=init.RR;

% free up memory
%---------------
clear data_info syst init

is_endogenous_switching=false;
if iscell(Q)
    Q0=Q{1};
    % check the possibility of endogenous probabilities
    is_endogenous_switching=numel(Q)>1 && ~isempty(Q{2});
    if is_endogenous_switching
        transition_matrix=Q{2};
        Vargs={};
        if numel(Q)>2
            Vargs=Q{3};
        end
    end
    Q=Q0;
end
PAI=transpose(Q)*PAItt;

% matrices' sizes
%----------------
[p0,smpl]=size(data);
smpl=min(smpl,find(include_in_likelihood,1,'last'));
m=size(T{1},1);
h=numel(T);
nshocks=size(R{1},2);
c_last=0;if ~isempty(state_trend),c_last=size(state_trend{1},2);end
d_last=0;if ~isempty(data_trend),d_last=size(data_trend{1},2);end
rqr_last=size(RR{1},3);
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

% free up memory
%---------------
% M=rmfield(M,{'data','T','data_structure','include_in_likelihood',...
%     'data_trend','state_trend','t_dc_last','obs_id','a','P','Q','PAI00','H','RR'});

% definitions and options
%------------------------
twopi=2*pi;
nsteps=options.kf_nsteps;
riccati_tol=options.kf_riccati_tol;
kalman_tol=options.kf_tol;
store_filters=options.kf_filtering_level;

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
    occur=data_structure(:,t);
    p=sum(occur); % number of observables to be used in likelihood computation
    y=data(occur,t);
    obsOccur=obs_id(occur); %<-- Z=M.Z(occur,:);
    
    likt=0;
    for st=1:h
        % forecast of observables
        %------------------------
        yf=a{st}(obsOccur); %<-- yf=Z*a{st};
        
        if d_last>0 % <-- ~isempty(data_trend)
            yf=yf+data_trend{st}(occur,min(t,d_last));
        end
        
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
            P{st}=P{st}-K(:,occur,st)*PZt';
            
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
    
    if store_filters>1
        store_updates();
    end
    
    % Likelihood computation
    %-----------------------
    Incr(t)=log(likt);
    
    % endogenous probabilities (conditional on time t information)
    %-------------------------------------------------------------
    att=a;
    if ~is_steady
        Ptt=P;
    end
    if is_endogenous_switching
        [Q,retcode]=online_function_evaluator(transition_matrix,att{1},Vargs{:});
        if retcode
            return
        end
    end
    
    % Probabilities predictions
    %--------------------------
    if h>1
        PAI=Q'*PAItt;
    end
    
    % state and state covariance prediction
    %--------------------------------------
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
        a{splus}=T{splus}(:,any_T)*a{splus}(any_T); % a{splus}=T{splus}*a{splus};
        if ~is_steady
            P{splus}=T{splus}(:,any_T)*P{splus}(any_T,any_T)*Tt{splus}(any_T,:)+RR{splus}(:,:,min(t,rqr_last));
            P{splus}=symmetrize(P{splus});
%             P{splus}=T{splus}*P{splus}*Tt{splus}+RR{splus}(:,:,min(t,rqr_last));
        end
        
        if c_last>0 % <-- ~isempty(state_trend)
            a{splus}=a{splus}+state_trend{splus}(:,min(t+1,c_last));
        end
    end
    
    if store_filters>0
        store_predictions()
    end
    
    if ~is_endogenous_switching && ~is_steady % && h==1
        is_steady=check_steady_state_kalman(is_steady);
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
                smoothing_step(Filters.a{s0}(:,1,t),r(:,s0),...
                K_store{s0}(:,occur,t),Filters.P{s0}(:,:,t),T{s0},R{s0},Z,...
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


    function is_steady=check_steady_state_kalman(is_steady)
        % the likelihood affects the values of the updated
        % probabilities. In turn, the updated probabilities affect the
        % values of the predicted probabilities. The predicted
        % probabilities enter the collapsing of the covariances and so,
        % we never reach the steady state in this case. So far I am
        % 100% sure about this. But in order to be 101% sure, I would
        % like to run an example.
        if t>no_more_missing
            discrep=max(abs(K(:)-oldK));
            is_steady=discrep<riccati_tol;
%             fprintf(1,'iteration %4.0f  discrepancy %4.4f\n',t,discrep);
        end
        oldK=K(:);
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
                Filters.a{splus_}(:,istep_,t+1)=T{splus_}*Filters.a{splus_}(:,istep_-1,t+1);
                if c_last>0 % <-- ~isempty(state_trend)
                    Filters.a{splus_}(:,istep_,t+1)=Filters.a{splus_}(:,istep_,t+1)+state_trend{splus_}(:,min(t+istep_,c_last));
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
                Filters.a{state}(:,1,1)=a{state};
                Filters.P{state}(:,:,1)=P{state};
            end
            Filters.PAI=zeros(h,smpl+1);
            Filters.PAI(:,1)=PAI;
            for istep=2:nsteps
                for state=1:h
                    Filters.a{state}(:,istep,1)=T{state}*Filters.a{state}(:,istep-1,1);
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

% Symmetrizer
%------------
function A=symmetrize(A)
    A=.5*(A+A.');
end

% Smoother
%------------
function [atT,etat,rlag]=smoothing_step(a,r,K,P,T,R,Z,iF,v)
L=T-T*K*Z;
% Note that Durbin and Koopman define K=T*P*Z'*iF, while here it is defined
% as K=P*Z'*iF. Hence, in the definition of Lt, I have to premultiply K by
% T
rlag=Z'*iF*v+L'*r;
atT=a+P*rlag;
% Q=eye(exo_nbr) in this algorithm...
etat=R'*rlag; % <--- etat=Rt'*rt;
% The state equation in Durbin and Koopman is
% a_{t+1}=T_{t}*a_{t}+R_{t}*eta_{t}, whereas the state equation in our
% models is a_{t}=T_{t}*a_{t-1}+R_{t}*eta_{t}, and this explains the change
% I made to the shock smoothing equation. With this change, if we define an
% auxiliary variable to be equal to a shock, we should retrieve the same
% values for both in the smoothing. Without the change, there will be a
% "delay"
end