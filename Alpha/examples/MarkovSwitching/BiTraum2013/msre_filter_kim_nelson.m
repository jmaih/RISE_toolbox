function [loglik,Incr,retcode,Filters]=msre_filter_kim_nelson(M,options)

% this filter assumes a state space of the form
% X_t=c_t(:,t,st)+T(:,:,st)*X_{t-1}+R(:,:,st)*eps_t
% y_t=d_t(:,t,st)+Z*X_t+eta_t
% where c_t(:,t,st) and d_t(:,t,st) are, possibly time-varying, deterministic terms
% the covariance matrices can be time-varying

% Symmetrizer
%------------
symmetrize=M.symmetrize;

% add fields if necessary
%------------------------
M=add_fields(M);

% data
%-----
data=M.data;
data_structure=M.data_structure;
include_in_likelihood=M.include_in_likelihood;
no_more_missing=M.no_more_missing;

% state matrices
%---------------
T=M.T;
R=M.R;
RR=M.RR;
H=M.H;
d=M.d;
c=M.c;
obs_id=M.obs_id;

% initial conditions
%-------------------
a=M.a;
P=M.P;
Q=M.Q0;
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
PAItt=M.PAI00;
PAI=transpose(Q)*PAItt;

% matrices' sizes
%----------------
[p0,smpl]=size(data);
[m,~,h]=size(T);
nshocks=size(R,2);
c_last=0;if ~isempty(c),c_last=size(c,2);end
d_last=0;if ~isempty(d),d_last=size(d,2);end
rqr_last=size(RR,3);
h_last=0;
if ~isempty(H)
    h_last=size(H,3);
end

% few transformations
%--------------------
any_T=false(1,m);
for st=1:h
    any_T=any_T|any(abs(T(:,:,st))>1e-9,1);
end
Tt=permute(T,[2,1,3]); % transpose

% free up memory
%---------------
% M=rmfield(M,{'data','T','data_structure','include_in_likelihood',...
%     'd','c','t_dc_last','obs_id','a','P','Q','PAI00','H','RR'});

% definitions and options
%------------------------
twopi=2*pi;
nsteps=options.nsteps;
riccati_tol=options.riccati_tol;
kalman_tol=options.kalman_tol;
store_filters=options.store_filters;

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
iF=zeros(p0,p0,h);
K=zeros(m,p0,h);
v=zeros(p0,h);

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
    p=sum(occur);
    y=data(occur,t);
    obsOccur=obs_id(occur); %<-- Z=M.Z(occur,:);
    
    likt=0;
    for st=1:h
        % forecast of observables
        %------------------------
        yf=a(obsOccur,st); %<-- yf=Z*a{st};
        
        if ~isempty(d)
            yf=yf+d(occur,min(t,d_last),st);
        end
        
        % forecast errors and variance
        %-----------------------------
        v(occur,st)=y-yf;
        
        if ~is_steady
            PZt=P(:,obsOccur,st); % PZt=<-- P{st}*Z';
            
            Fst=PZt(obsOccur,:); % <--- F=Z*PZt+H{st}(occur,occur);
            if h_last>0
                Fst=Fst+H(occur,occur,min(t,h_last),st);
            end
            
            % Kalman gain (for update)
            %-------------------------
            % number of observables to be used in likelihood computation
            detF=det(Fst);
            failed=detF<=0;
            if ~failed
                iF(occur,occur,st)=Fst\eye(p);
                failed=any(any(isnan(iF(occur,occur,st))));
            end
            if failed
                retcode=305;
                return
            end
            K(:,occur,st)=PZt*iF(occur,occur,st); % K=PZt/F{st};
            
            % state covariance update (Ptt=P-P*Z*inv(F)*Z'*P)
            %------------------------------------------------
            P(:,:,st)=P(:,:,st)-K(:,occur,st)*PZt';
            
            twopi_p_dF(st)=twopi^p*det(Fst);
        end
        % state update (att=a+K*v)
        %-------------------------
        a(:,st)=a(:,st)+K(:,occur,st)*v(occur,st);
        
        f01=(twopi_p_dF(st)*exp(v(occur,st)'*iF(occur,occur,st)*v(occur,st)))^(-0.5);
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
        [Q,retcode]=online_function_evaluator(transition_matrix,att(:,1),Vargs{:});
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
    if ~is_steady
        P=zeros(m,m,h);
    end
    a=zeros(m,h);
    for splus=1:h
        for st=1:h
            if h==1
                pai_st_splus=1;
            else
                pai_st_splus=Q(st,splus)*PAItt(st)/PAI(splus);
            end
            a(:,splus)=a(:,splus)+pai_st_splus*att(:,st);
            if ~is_steady
                P(:,:,splus)=P(:,:,splus)+pai_st_splus*Ptt(:,:,splus);
            end
        end
        a(:,splus)=T(:,any_T,splus)*a(any_T,splus);
        if ~is_steady
            P(:,:,splus)=T(:,any_T,splus)*P(any_T,any_T,splus)*Tt(any_T,:,splus)+RR(:,:,min(t,rqr_last),splus);
            P(:,:,splus)=symmetrize(P(:,:,splus));
%             P(:,:,splus)=T(:,:,splus)*P(:,:,splus)*Tt(:,:,splus)+RR(:,:,min(t,rqr_last),splus);
        end
        
        if ~isempty(c)
            a(:,splus)=a(:,splus)+c{splus}(:,min(t+1,c_last));
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
    for t=smpl:-1:1
        Q=Filters.Q(:,:,t);
        occur=data_structure(:,t);
        Z=M.Z(occur,:);
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
                K_store{s0}(:,occur,t),Filters.P{s0}(:,:,t),T(:,:,s0),...
                R(:,:,s0),Z,...
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
            Filters.att{st_}(:,1,t)=a(:,st_);
            Filters.Ptt{st_}(:,:,t)=P(:,:,st_);
            if store_filters>2
                K_store{st_}(:,occur,t)=K(:,occur,st_);
                iF_store{st_}(occur,occur,t)=iF(occur,occur,st_);
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
                if ~isempty(c)
                    Filters.a{splus_}=Filters.a{splus_}+c(:,min(t+istep_,c_last),splus_);
                end
            end
        end
    end
    function Filters=initialize_storage()
        Filters=struct();
        if store_filters>0 % store filters
            Filters.a=repmat({zeros(m,nsteps,smpl+1)},[h,h]);
            Filters.P=repmat({zeros(m,m,smpl+1)},[h,h]);
            for state=1:h
                for state2=1:h
                    Filters.a{state,state2}(:,1,1)=a(:,state,state2);
                    Filters.P{state,state2}(:,:,1)=P(:,:,state,state2);
                end
            end
            Filters.PAI=zeros(h,smpl+1);
            Filters.PAI(:,1)=PAI;
            for istep=2:nsteps
                for state=1:h
                    for state2=1:h
                        Filters.a{state,state2}(:,istep,1)=T(:,:,state)*Filters.a{state,state2}(:,istep-1,1);
                    end
                end
            end
            Filters.Q=zeros(h,h,smpl+1);
            Filters.Q(:,:,1)=Q;
            if store_filters>1 % store updates
                Filters.att=repmat({zeros(m,1,smpl)},h,h);
                Filters.Ptt=repmat({zeros(m,m,smpl)},h,h);
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

function M=add_fields(M)
myfields={
    'c',[]
    'd',[]
    'Q0',1 % initial transition matrix
    'PAI00',1 % initial state probabilities
    'H',[] % covariance matrix of measurement errors
    };
for jj=1:size(myfields,1)
    field=myfields{jj,1};
    if ~isfield(M,field)
        M.(field)=myfields{jj,2};
    end
end

end