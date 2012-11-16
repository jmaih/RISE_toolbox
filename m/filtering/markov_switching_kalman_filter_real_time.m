function [LogLik,Incr,retcode,Filters] = markov_switching_kalman_filter_real_time(obs_id,...
    y,... % data
    rest_id,... % restrictions location
    MU,... % conditions to hit
    OMG,... % covariance matrix (uncertainty) of the conditions
    Hypothesis,... % how to forecast
    T,...
    R,...
    steady_state,...
    Q,... % transition matrix
    H,... % covariance of measurement errors
    Options,...
    WB)

% all rows of Q should sum to 1

defaults=markov_switching_kalman_filter();
if nargin==0
    if nargout>1
        error([mfilename,':: with no input argument, the number of output arguments cannot exceed 1'])
    end
    LogLik=defaults;
    return
end
% all rows of Q should sum to 1
kf_algorithm='';
kf_tol=0;
kf_filtering_level=0;

try
    narginchk(11,13)
catch %#ok<CTCH>
    % for backward compatibility
    error(nargchk(11,13,nargin,'struct')) %#ok<NCHKN>
end

if nargin<13
    WB=[];
    if nargin<12
        Options=[];
        if nargin<10
            H=[];
        end
    end
end

data_occurrence=~isnan(y);

filtering_fields=fieldnames(defaults);
for ii=1:numel(filtering_fields)
    v=filtering_fields{ii};
    if isfield(Options,v)
        defaults.(v)=Options.(v);
    end
    eval([v,'=defaults.(v);'])
end
init_options=kalman_initialization();
init_fields=fieldnames(init_options);
for ii=1:numel(init_fields)
    v=init_fields{ii};
    init_options.(v)=defaults.(v);
end


deterministic=false;
if ~isempty(WB)
    W=WB{1}; %R_det
    B=WB{2}; %exo_data
    clear WB
    deterministic=true;
end

endogenous_switching=~isempty(Q{2});
Q0=Q{1};
if endogenous_switching
    transition_matrix=Q{2};
    Vargs=Q{3};
end
clear Q


ExpandedFlag=true;
% initialize this in case the kf_filtering_level does not go to the end
Filters=[];
Incr=[];
LogLik=nan;

[endo_nbr,~,order,h]=size(R);
% currently this would not work for nsteps>1 but we will see.

[pp,smpl]=size(y);
if isempty(H)
    H=zeros(pp,pp,h);
end

if isempty(MU)
    MU=nan(0,order,smpl);
end
[nrest,number_of_conditioning_periods,smpl0]=size(MU);
% Build shock impact
if ~isequal(numel(rest_id),nrest)
    error([mfilename,':: number of restrictions inconsistent with size of MU'])
end
if ~isequal(smpl0,smpl)
    error([mfilename,':: sample lengths in both y and MU should be the same when MU is not empty'])
end

% expansion order has to be consistent with the number of steps taken...
Nsteps=min(order,number_of_conditioning_periods);
DPHI=cell(1,h);
DT=cell(1,h);
for st=1:h
    [DPHI{st},DT{st}]=BuildShockRestrictions(...
        T(:,:,st),R(:,:,:,st),rest_id,[],number_of_conditioning_periods,Nsteps,Hypothesis);
end
shock_span=size(DPHI{1},2);

mm=endo_nbr+shock_span;
steady_state=[steady_state;zeros(shock_span,h)];

switch lower(deblank(kf_algorithm))
    case 'lwz'
        hstar=1;
        filt_algo=1;
    case 'kn'
        hstar=h;
        filt_algo=2;
    otherwise
end

%% Initialization
a01=zeros(mm,h,hstar);
for s1=1:h
    for s0=1:hstar
        a01(1:endo_nbr,s1,s0)=steady_state(:,s1);
        if deterministic
            occur_x=~isnan(B(:,1));
            d1=W(:,occur_x,s1)*B(occur_x,1); % first observation of exogenous
            a01(1:endo_nbr,s1,s0)=a01(1:endo_nbr,s1,s0)+d1;
        end
    end
end
P01=zeros(mm,mm,h,hstar);
for s1=1:h
    [~,tmp,PAItt,start,retcode]=kalman_initialization(...
        T(:,:,s1),R(:,:,1,s1)*R(:,:,1,s1)',Q0,init_options);
    if retcode==0
        for s0=1:hstar
            P01(1:endo_nbr,1:endo_nbr,s1,s0)=tmp;
            P01(endo_nbr+1:end,endo_nbr+1:end,s1,s0)=eye(shock_span);
        end
    else
        return
    end
end
PAI=transpose(Q0)*PAItt;

%% pre-allocation
v=zeros(pp,h,hstar);
F=zeros(pp,pp,h,hstar);
iF=zeros(pp,pp,h,hstar);
PAI01y=zeros(h,hstar);
K=zeros(mm,pp,h,hstar);
Tt=zeros(mm,mm,h);
Rt=zeros(mm,shock_span,h);

if  kf_filtering_level
    Filters.BIGPAI=zeros(h,smpl+1);
    Filters.BIGPAI(:,1)=PAI;
    a01_store=zeros(mm,h,hstar,smpl+1);
    a_store=zeros(mm,h,smpl+1);
    P_store=zeros(mm,mm,h,hstar,smpl+1);
    for s1=1:h
        for s0=1:hstar
            P_store(:,:,s1,s0,1)=P01(:,:,s1,s0);
        end
    end
    if  kf_filtering_level>1
        Filters.BIGPAI_tt=zeros(h,smpl);
        att_store=zeros(mm,h,smpl);
        if  kf_filtering_level>2
            v_store=zeros(pp,h,hstar,smpl);
            iF_store=zeros(pp,pp,h,hstar,smpl);
            K_store=zeros(mm,pp,h,hstar,smpl);
            T_store=zeros(mm,mm,h,smpl);
            R_store=zeros(mm,shock_span,h,smpl);
        end
    end
    Q=nan(h,h,smpl);
end


%% recursions
Record=cell(1,hstar);
Incr=nan(smpl,1);
Qt=Q0;
for t=1:smpl
    likt=0;
    occur=data_occurrence(:,t);
    for s1=h:-1:1
        data=y(:,t);
        % the substracting of the steady state really makes more sense when
        % the there is only one regime or whether the steady state is the
        % same in all the regimes...
        for s0=hstar:-1:1
            %% forecast errors
            v(occur,s1,s0)=data(occur)-a01(obs_id(occur),s1,s0);
            % Symmetrize it first, just in case
            F(occur,occur,s1,s0)=P01(obs_id(occur),obs_id,s1(occur),s0)+H(occur,occur,s1);% symmetrize()
            [ispd,dF,iF10,F(occur,occur,s1,s0)]=CheckPositiveDefiniteness(F(occur,occur,s1,s0));
            if ~ispd
                retcode=305;
                return
            end
            iF(occur,occur,s1,s0)=iF10;
            f01=conditional_likelihood(v(occur,s1,s0),iF(occur,occur,s1,s0),dF,sum(occur));
            if filt_algo==1
                pai01=PAI(s1);
            elseif filt_algo==2
                pai01=Qt(s0,s1)*PAItt(s0);
            end
            PAI01y(s1,s0)=pai01*f01;
            likt=likt+PAI01y(s1,s0);
        end
    end
    
    PAI01_tt=PAI01y/likt;
    if likt<kf_tol && (any(isnan(PAI01_tt))||any(isinf(PAI01_tt)))
        retcode=306;
        return
    end
    PAItt=sum(PAI01_tt,2);
    Incr(t)=log(likt);
    
    %% Updating
    if filt_algo==1
        att=zeros(mm,h);
        Ptt=zeros(mm,mm,h);
        for s1=1:h
            [att(:,s1),Ptt(:,:,s1),K(:,occur,s1)]=kalman_update(a01(:,s1),P01(:,:,s1),iF(occur,occur,s1),v(occur,s1),obs_id(occur));
        end
    elseif filt_algo==2
        % update and collapse
        [att,Ptt,K(:,occur,:,:)]=update_and_collapse(a01,P01,iF,v,obs_id,mm,pp,h,...
            hstar,PAI01_tt,PAItt);
    end
    
    %% prediction
    PAI=Qt'*PAItt;
    OMGt=OMG;
    dt_plus_1=0;
    for s1=h:-1:1
        MUt=MU(:,:,t);
        if deterministic && t<smpl
            % the last prediction will be strongly biased since we don't
            % have out-of-sample data
            occur_x=~isnan(B(:,t+1));
            dt_plus_1=W(:,occur_x,s1)*B(occur_x,t+1);
            if ~isempty(MUt)
                % change the first period only: substract a term that will
                % be added back
                MUt(:,1)=MUt(:,1)-dt_plus_1(rest_id);
            end
        end
        for s0=hstar:-1:1
            if filt_algo==1
                [a01(:,s1),P01(:,:,s1),Tt(:,:,s1),Rt(:,:,s1),Record{s1}]=predict_while_collapsing(...
                    att,Ptt,Qt,PAI,PAItt,...
                    MUt,OMGt,DPHI{s1},DT{s1},Record{s1},ExpandedFlag,...
                    s1);
                a01(:,s1)=a01(:,s1)+dt_plus_1;
            elseif filt_algo==2
                [a01(:,s1,s0),P01(:,:,s1,s0),Tt(:,:,s1),Rt(:,:,s1),Record{s1}]=kalman_prediction(...
                    T(:,:,s1),R(:,:,:,s1),att(:,s0),Ptt(:,:,s0),...
                    MUt,OMGt,DPHI{s1},DT{s1},Record{s1},ExpandedFlag);
                a01(:,s1,s0)=a01(:,s1,s0)+dt_plus_1;
            end
        end
    end
    
    if endogenous_switching
        [Qtmp,retcode]=kron(transition_matrix,att(:,1),Vargs{:});
%         [Qtmp,retcode]=transition_matrix(att(:,1),Vargs{:});
        if retcode
            return
        end
    else
        Qtmp=Q0;
    end
    
    %% store for kf_filtering_level
    if kf_filtering_level
        Filters.BIGPAI(:,t+1)=PAI;
        P_store(:,:,:,:,t+1)=P01;
        a01_store(:,:,:,t+1)=a01;
        for s1=1:h
            for s0=1:hstar
                pai01_joint=Qt(s0,s1)*PAItt(s1);
                a_store(:,s1,t+1)=a_store(:,s1,t+1)+pai01_joint*a01(:,s1,s0);
            end
            a_store(:,s1,t+1)=a_store(:,s1,t+1)/PAItt(s1);
        end
        if  kf_filtering_level>1
            Filters.BIGPAI_tt(:,t)=PAItt;
            att_store(:,:,t)=att;
            if  kf_filtering_level>2
                v_store(:,:,:,t)=v;
                iF_store(:,:,:,:,t)=iF;
                K_store(:,:,:,:,t)=K;
                T_store(:,:,:,t)=Tt;
                R_store(:,:,:,t)=Rt;
            end
        end
        Q(:,:,t)=Qtmp;
    end
    Qt=Qtmp;
end
LogLik=sum(Incr(start:end));

if kf_filtering_level
    if kf_filtering_level>1
        if kf_filtering_level>2
            r=zeros(mm,h,hstar);
            pai_lead_now=zeros(h);
            Zt=eye(mm);
            Zt=Zt(obs_id,:);
            Filters.BIGPAI_tT=zeros(h,smpl);
            for t=smpl:-1:1
                Qt=Q(:,:,t);
                occur=data_occurrence(:,t);
                for s0=1:h
                    for slead=1:h
                        % joint probability of being in state slead tomorrow and
                        % beeing in state s0 today, conditional on all the
                        % information available
                        if t==smpl
                            pai_lead_now(s0,slead)=Qt(s0,slead)*Filters.BIGPAI_tt(s0,t);
                        else
                            pai_lead_now(s0,slead)=Qt(s0,slead)*Filters.BIGPAI_tt(s0,t)*Filters.BIGPAI_tT(slead,t+1)/Filters.BIGPAI(slead,t+1);
                        end
                        % smoothed probabilities
                        Filters.BIGPAI_tT(s0,t)=Filters.BIGPAI_tT(s0,t)+pai_lead_now(s0,slead);
                    end
                    % smoothed state
                    %             [alphat,etat,rlag]=smoothing_step(at,rt,Kt,Pt,Tt,Rt,Zt,iFt,vt)
                    if filt_algo==1
                        [Filters.alphat(:,s0,t),Filters.eta(:,s0,t),r(:,s0)]=...
                            smoothing_step(a01_store(:,s0,1,t),r(:,s0),...
                            K_store(:,occur,s0,1,t),P_store(:,:,s0,1,t),T_store(:,:,s0,t),R_store(:,:,s0,t),Zt(occur,:),...
                            iF_store(occur,occur,s0,1,t),v_store(occur,s0,1,t));
                    elseif filt_algo==2
                        eta0=0;
                        alpha0=0;
                        for s1=1:hstar
                            [alphat01,eta01,r(:,s0,s1)]=smoothing_step(a01_store(:,s0,s1,t),r(:,s0,s1),...
                                K_store(:,occur,s0,s1,t),P_store(:,:,s0,s1,t),T_store(:,:,s1,t),R_store(:,:,s1,t),Zt(occur,:),...
                                iF_store(occur,occur,s0,s1,t),v_store(occur,s0,s1,t));
                            alpha0=alpha0+pai_lead_now(s0,s1)*alphat01;
                            eta0=eta0+pai_lead_now(s0,s1)*eta01;
                        end
                        Filters.alphat(:,s0,t)=alpha0/Filters.BIGPAI_tT(s0,t);
                        Filters.eta(:,s0,t)=eta0/Filters.BIGPAI_tT(s0,t);
                    end
                    % measurement errors
                    Filters.epsilon(1:pp,s0,t)=0;
                    Filters.epsilon(occur,s0,t)=y(occur,t)-Filters.alphat(obs_id(occur),s0,t);
                end
                % correction for the smoothed probabilities [the approximation involved does not always work
                % especially when dealing with endogenous switching.
                SumProbs=sum(Filters.BIGPAI_tT(:,t));
                if abs(SumProbs-1)>1e-8
                    Filters.BIGPAI_tT(:,t)=Filters.BIGPAI_tT(:,t)/SumProbs;
                    if filt_algo==2
                        Filters.alphat(:,:,t)=Filters.alphat(:,:,t)*SumProbs;
                        Filters.eta(:,:,t)=Filters.eta(:,:,t)*SumProbs;
                    end
                end
            end
            % now in real-time, the smoothed shocks of interest are to be found
            % after the endogenous variables in the state vector and not in the
            % usual smoothed disturbances. This is because in the presence of
            % conditional information, smoothed shocks are no longer IID whereas
            % smoothed disturbances are. In order to be able to reconstruct the
            % smoothed series in the absence of the time varying state matrices, we
            % need the smoothed shocks, not the smoothed disturbances. And so we do
            % the following
            Filters.eta=Filters.alphat(endo_nbr+1:end,:,:);
            % now we can resize the smoothed variables array
            Filters.alphat=Filters.alphat(1:endo_nbr,:,:);
            % and in fact we need to resize all the other filters as well since we
            % don't need the filtered shocks and the updated shocks
        end % if kf_filtering_level>2
        Filters.att=att_store;
        Filters.att=Filters.att(1:endo_nbr,:,:);
        
    end % if kf_filtering_level>1
    
    Filters.a=a_store;
    Filters.a=Filters.a(1:endo_nbr,:,:);
    % It can be checked that in the absence of conditional information,
    % smoothed disturbances and smoothed shocks are identical.
end % if kf_filtering_level

%=========================== Sub-function ===========================================
    function [a,P,Tt,Rt,Record]=predict_while_collapsing(att,Ptt,Qt,PAI,PAItt,...
            MUt,OMGt,DPHI,DT,Record,ExpandedFlag,slead)
        a=0;
        P=0;
        Tt=0;
        Rt=0;
        for snow=1:h
            pai_snow_Over_slead=Qt(snow,slead)*PAItt(snow)/PAI(slead);
            [a00_,P00_,Tt_snow,Rt_snow,Record]=kalman_prediction(...
                T(:,:,slead),R(:,:,:,slead),att(:,snow),Ptt(:,:,snow),...
                MUt,OMGt,DPHI,DT,Record,ExpandedFlag);
            a=a+pai_snow_Over_slead*a00_;
            P=P+pai_snow_Over_slead*P00_;
            Tt=Tt+pai_snow_Over_slead*Tt_snow;
            Rt=Rt+pai_snow_Over_slead*Rt_snow;
        end
        %%%        P=symmetrize(P);
    end
end

