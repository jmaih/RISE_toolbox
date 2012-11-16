function [LogLik,Incr,retcode,Filters] = markov_switching_kalman_filter(obs_id,y,T,...
    R,steady_state,Q,H,Options,WB)
%  Detailed explanation to come here

% all rows of Q should sum to 1
kf_algorithm='';
kf_tol=0;
kf_filtering_level=0;
filt_options=struct('kf_algorithm','lwz',...%     alternative is kn (Kim and Nelson)
    'kf_tol',1.0000e-020,...
    'kf_filtering_level',3);
% 0: no filters,
% 1: filtering,
% 2: filtering+updating,
% 3: filtering+updating+smoothing
% this decides the initial condition for both the state and the markov chain distributions in the
% the kalman filter.
% now add the defaults for the initialization process
init_options=kalman_initialization();
defaults=mergestructures(filt_options,init_options);
if nargin==0
    if nargout>1
        error([mfilename,':: with no input argument, the number of output arguments cannot exceed 1'])
    end
    LogLik=defaults;
    return
end

try
    narginchk(6,9)
catch %#ok<CTCH>
    % for backward compatibility
    error(nargchk(6,9,nargin)) %#ok<NCHKN>
end

if nargin<9
    WB=[];
    if nargin<8
        Options=[];
        if nargin<7
            H=[];
        end
    end
end
data_occurrence=~isnan(y);

filtering_fields=fieldnames(filt_options);
for ii=1:numel(filtering_fields)
    v=filtering_fields{ii};
    if isfield(Options,v)
        eval([v,'=Options.(v);'])
    else
        eval([v,'=defaults.(v);'])
    end
end
init_fields=fieldnames(init_options);
for ii=1:numel(init_fields)
    v=init_fields{ii};
    if isfield(Options,v)
        init_options.(v)=Options.(v);
    end
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

% initialize this in case the filtering does not go to the end
Filters=[];
Incr=[];
LogLik=nan;

[endo_nbr,~,nsteps,h]=size(R);

if nsteps>1
    error([mfilename,':: this seems to be a candidate for the real-time filter'])
end

[pp,smpl]=size(y);
if isempty(H)
    H=zeros(pp,pp,h);
end

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
a01=zeros(endo_nbr,h,hstar);
for s1=1:h
    for s0=1:hstar
        a01(:,s1,s0)=steady_state(:,s1);
        if deterministic
            occur_x=~isnan(B(:,1));
            d1=W(:,occur_x,s1)*B(occur_x,1); % first observation of exogenous
            a01(:,s1,s0)=a01(:,s1,s0)+d1; %
        end
    end
end
P01=zeros(endo_nbr,endo_nbr,h,hstar);
RR=zeros(endo_nbr,endo_nbr,h);

for s1=1:h
    RR(:,:,s1)=R(:,:,s1)*transpose(R(:,:,s1));
    % the correct form is RR(:,:,s1)=R(:,:,:,s1)*transpose(R(:,:,:,s1));
    % but still the thing above works since the 3rd dimension is unity
    
    % initialize the filter but also modify T and RR if there are
    % nonstationary variables
    [~,tmp,PAItt,start,retcode]=kalman_initialization(T(:,:,s1),RR(:,:,s1),...
        Q0,init_options);
    if retcode==0
        for s0=1:hstar
            P01(:,:,s1,s0)=tmp;
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
K=zeros(endo_nbr,pp,h,hstar);

if  kf_filtering_level
    Filters.BIGPAI=zeros(h,smpl+1);
    Filters.BIGPAI(:,1)=PAI;
    a01_store=zeros(endo_nbr,h,hstar,smpl+1);
    a_store=zeros(endo_nbr,h,smpl+1);
    P_store=zeros(endo_nbr,endo_nbr,h,hstar,smpl+1);
    for s1=1:h
        for s0=1:hstar
            P_store(:,:,s1,s0,1)=P01(:,:,s1,s0);
        end
    end
    if  kf_filtering_level>1
        Filters.BIGPAI_tt=zeros(h,smpl);
        att_store=zeros(endo_nbr,h,smpl);
        if  kf_filtering_level>2
            v_store=zeros(pp,h,hstar,smpl);
            iF_store=zeros(pp,pp,h,hstar,smpl);
            K_store=zeros(endo_nbr,pp,h,hstar,smpl);
        end
    end
    Q=nan(h,h,smpl);
end

Incr=nan(smpl,1);
%% recursions
Qt=Q0;
for t=1:smpl
    likt=0;
    occur=data_occurrence(:,t);
    for s1=1:h
        data=y(:,t);
        % the substracting of the steady state really makes more sense when
        % the there is only one regime or whether the steady state is the
        % same in all the regimes...
        for s0=1:hstar
            %% forecast errors
            v(occur,s1,s0)=data(occur,1)-a01(obs_id(occur),s1,s0);
            % Symmetrize it first, just in case
            F(occur,occur,s1,s0)=P01(obs_id(occur),obs_id(occur),s1,s0)+H(occur,occur,s1); %symmetrize()
            [ispd,dF,iF10]=CheckPositiveDefiniteness(F(occur,occur,s1,s0));
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
        att=zeros(endo_nbr,h);
        Ptt=zeros(endo_nbr,endo_nbr,h);
        for s1=1:h
            [att(:,s1),Ptt(:,:,s1),K(:,occur,s1)]=kalman_update(a01(:,s1),P01(:,:,s1),iF(occur,occur,s1),v(occur,s1),obs_id(occur));
        end
    elseif filt_algo==2
        % update and collapse
        [att,Ptt,K(:,occur,:,:)]=update_and_collapse(a01,P01,iF(occur,occur,:,:),...
            v(occur,:,:),obs_id(occur),endo_nbr,sum(occur),h,...
            hstar,PAI01_tt,PAItt);
    end
    
    %% prediction
    PAI=Qt'*PAItt;
    dt_plus_1=0;
    for s1=1:h
        if deterministic && t<smpl
            % the last prediction will be strongly biased since we don't
            % have out-of-sample data
            occur_x=~isnan(B(:,t+1));
            dt_plus_1=W(:,occur_x,s1)*B(occur_x,t+1);
        end
        for s0=1:hstar
            if filt_algo==1
                [a01(:,s1),P01(:,:,s1)]=predict_while_collapsing(att,Ptt,Qt,PAI,PAItt,s1);
                a01(:,s1)=a01(:,s1)+dt_plus_1;
            elseif filt_algo==2
                [a01(:,s1,s0),P01(:,:,s1,s0)]=kalman_prediction(T(:,:,s1),RR(:,:,s1),...
                    att(:,s0)-steady_state(:,s0),Ptt(:,:,s0));
                % PAI(s1)=PAI(s1)+Qt(s0,s1)*PAItt(s0);
                a01(:,s1,s0)=a01(:,s1,s0)+steady_state(:,s1)+dt_plus_1;
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
        a01_store(:,:,:,t+1)=a01;
        for s1=1:h
            for s0=1:hstar
                pai01_joint=Qt(s0,s1)*PAItt(s1);
                a_store(:,s1,t+1)=a_store(:,s1,t+1)+pai01_joint*a01(:,s1,s0);
            end
            a_store(:,s1,t+1)=a_store(:,s1,t+1)/PAItt(s1);
        end
        Filters.BIGPAI(:,t+1)=PAI;
        P_store(:,:,:,:,t+1)=P01;
        if  kf_filtering_level>1
            Filters.BIGPAI_tt(:,t)=PAItt;
            att_store(:,:,t)=att;
            if  kf_filtering_level>2
                v_store(:,:,:,t)=v;
                iF_store(:,:,:,:,t)=iF;
                K_store(:,:,:,:,t)=K;
            end
        end
        Q(:,:,t)=Qtmp;
    end
    Qt=Qtmp;
end
LogLik=sum(Incr(start:end));

if kf_filtering_level
    Filters.a=a_store;
    if kf_filtering_level>1
        Filters.att=att_store;
        
        if kf_filtering_level>2
            r=zeros(endo_nbr,h,hstar);
            pai_lead_now=zeros(h);
            Zt=eye(endo_nbr);
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
                    % [alphat,etat,rlag]=smoothing_step(at,rt,Kt,Pt,Tt,Rt,Zt,iFt,vt)
                    if filt_algo==1
                        [Filters.alphat(:,s0,t),Filters.eta(:,s0,t),r(:,s0)]=...
                            smoothing_step(a01_store(:,s0,1,t),r(:,s0),...
                            K_store(:,occur,s0,1,t),P_store(:,:,s0,1,t),T(:,:,s0),R(:,:,s0),Zt(occur,:),...
                            iF_store(occur,occur,s0,1,t),v_store(occur,s0,1,t));
                    elseif filt_algo==2
                        eta0=0;
                        alpha0=0;
                        for s1=1:hstar
                            [alphat01,eta01,r(:,s0,s1)]=smoothing_step(a01_store(:,s0,s1,t),r(:,s0,s1),...
                                K_store(:,occur,s0,s1,t),P_store(:,:,s0,s1,t),T(:,:,s1),R(:,:,s1),Zt(occur,:),...
                                iF_store(occur,occur,s0,s1,t),v_store(occur,s0,s1,t));
                            alpha0=alpha0+pai_lead_now(s0,s1)*alphat01;
                            eta0=eta0+pai_lead_now(s0,s1)*eta01;
                        end
                        Filters.alphat(:,s0,t)=alpha0/Filters.BIGPAI_tT(s0,t);
                        Filters.eta(:,s0,t)=eta0/Filters.BIGPAI_tT(s0,t);
                    end
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
        end % if kf_filtering_level>2
    end % if kf_filtering_level>1
end % if kf_filtering_level

%=========================== Sub-function ===========================================
    function [a,P]=predict_while_collapsing(att,Ptt,Qt,PAI,PAItt,slead)
        a=0;
        P=0;
        pai_snow_Over_slead=1;
        for snow=1:h
            if h>1
                pai_snow_Over_slead=Qt(snow,slead)*PAItt(snow)/PAI(slead);
            end
            [a00,P00]=kalman_prediction(T(:,:,slead),RR(:,:,slead),...
                att(:,snow)-steady_state(:,snow),Ptt(:,:,snow));
            a=a+pai_snow_Over_slead*a00;
            P=P+pai_snow_Over_slead*P00;
        end
        a=a+steady_state(:,slead);
        %%%        P=symmetrize(P);
    end
end


