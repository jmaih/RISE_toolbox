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

% free up memory
%---------------
m=size(T{1},1);
h=numel(T);
npages=data_info.npages;
% nshocks=size(R{1},2);
kmax=size(R{1},3)-1;
horizon=kmax+1;
horizon=min(npages-1,horizon); % the first page is hard information
ncp=horizon;
hypothesis=options.forecast_conditional_hypothesis;

restr_y_id=data_info.restr_y_id;
restr_x_id=data_info.restr_x_id;
if ~isempty(restr_x_id)
    error('conditioning on shocks is not allowed in estimation')
end

DPHI=cell(1,h);
DT=cell(1,h);
for st=1:h
    [DPHI{st},DT{st}]=utils.forecast.conditional.build_shock_restrictions(...
        T{st},R{st},...
        restr_y_id,restr_x_id,...
        ncp,... data availability
        horizon,... horizon +or-1
        hypothesis);
    % redress: remove the third dimension
    %-------------------------------------
    R{st}=R{st}(:,:);
end

Q=Qfunc(a{1});
PAI=transpose(Q)*PAItt;

% matrices' sizes
%----------------
[p0,smpl]=size(data);
smpl=min(smpl,find(include_in_likelihood,1,'last'));
c_last=0;if ~isempty(state_trend),c_last=size(state_trend{1},2);end
% rqr_last=size(RR{1},3);
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
    MUt=data(restr_y_id,t,1+(1:horizon));
    
    % state and state covariance prediction
    %--------------------------------------
    for splus=1:h
        MUt_splus=MUt(:,:);
        if c_last>0 % remove the part that will be unexplained by the shocks
            MUt_splus=MUt_splus-state_trend{st}(restr_y_id,min(t+1,c_last))*ones(1,horizon);
        end
        [a{splus},P{splus},Tt{splus},Rt{splus},Record]=predict_while_collapsing(Record,expanded_flag);
        
        if c_last>0 % <-- ~isempty(state_trend)
            a{splus}(1:m)=a{splus}(1:m)+state_trend{splus}(:,min(t+1,c_last));
        end
    end
    
    if store_filters>0
        if store_filters>1
            store_updates();
        end
        store_predictions()
    end
    
    if ~is_steady % && h==1
        is_steady=check_steady_state_kalman(is_steady);
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
                smoothing_step(Filters.a{s0}(:,1,t),r(:,s0),...
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
        Filters.a{st}=Filters.a{st}(1:m,:,:); % second dimension is the real-time forecasting steps
        Filters.P{st}=Filters.P{st}(1:m,1:m,:);
        if store_filters>1
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

    function [a,P,Tt,Rt,Record]=predict_while_collapsing(Record,ExpandedFlag)
        a=0;
        P=0;
        Tt=0;
        Rt=0;
        for snow=1:h
            pai_snow_Over_slead=Q(snow,splus)*PAItt(snow)/PAI(splus);
            [a00_,P00_,Tt_snow,Rt_snow,Record]=kalman_prediction(...
                T{splus},R{splus},att{snow},Ptt{snow},...
                MUt_splus,OMGt,DPHI{splus},DT{splus},Record,ExpandedFlag);
            a=a+pai_snow_Over_slead*a00_;
            P=P+pai_snow_Over_slead*P00_;
            Tt=Tt+pai_snow_Over_slead*Tt_snow;
            Rt=Rt+pai_snow_Over_slead*Rt_snow;
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
            Filters.att{st_}(:,1,t)=att{st_};
            Filters.Ptt{st_}(:,:,t)=Ptt{st_};
            if expanded_flag
                K_store{st_}(:,occur,t)=K(:,occur,st_);
                iF_store{st_}(occur,occur,t)=iF{st_};
                v_store{st_}(occur,t)=v{st_};
                % state matrices
                Tt_store{st_}(:,:,t)=Tt{st_};
                Rt_store{st_}(:,:,t)=Rt{st_};
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
                Filters.a{splus_}(:,istep_,t+1)=Tt{splus_}*Filters.a{splus_}(:,istep_-1,t+1);
                if c_last>0 % <-- ~isempty(state_trend)
                    Filters.a{splus_}(1:m,istep_,t+1)=Filters.a{splus_}(1:m,istep_,t+1)+state_trend{splus_}(:,min(t+istep_,c_last));
                end
            end
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
end

function [a,P,Tt,Rt,Record]=kalman_prediction(T,R,att,Ptt,MUt,OMGt,DPHI,DT,Record,ExpandedFlag)
narginchk(4,10);
reduced=nargin==4;
bt=0;
Tt=T;
Rt=R;
Record=[];
if ~reduced
    reduced=reduced||(~ExpandedFlag && all(isnan(MUt(:))));
    if ~reduced
        [Tt,Rt,bt,~,Record]=utils.forecast.conditional.state_matrices(T,R,MUt,OMGt,DPHI,DT,Record,ExpandedFlag);
    end
end

RR=Rt*Rt';
% only compute the places where there is some action
test=true;
if test
    a=bt+Tt*att;
    P=Tt*Ptt*transpose(Tt)+RR;
else
    kk=any(Tt);
    a=bt+Tt(:,kk)*att(kk,:);
    P=Tt(:,kk)*Ptt(kk,kk)*transpose(Tt(:,kk))+RR;
end

% Make sure we remain symmetric
P=symmetrize(P);
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