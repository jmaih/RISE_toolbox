function [params,smoothed_vals]=gibbs_sampler(obj,params,smoothed_vals)
persistent y X priors pp n smpl Z initvals endo_ids stoch_vol_small...
    build_A build_OMG vectorize_OMG % vectorize_A
% initialize the gibbs
%---------------------
% output the smoothed values alongside the parameters...
if nargin==1
    params=struct();
    [y,X,n,smpl]=vartools.set_y_and_x(obj.data.y,obj.data.x,obj.nlags,obj.constant);
    if obj.time_varying_parameters(1)
        Z=build_Z();
    end
    
    % parameter priors
    %-----------------
    priors=obj.estimation.priors;
    params.x0=[priors.prior_mean]';
    
    % locations
    %----------
    %   endo_ids, exo_ids, params_ids, build_A, vectorize_A, buil_OMG, vectorize_OMG
    [endo_ids,~,pp,build_A,~,build_OMG,vectorize_OMG]=stochvol.format_blocks(obj);
    %     endo_ids=struct('y_id',y_id,'a_state_id',a_state_id,'a_matrix_id',a_matrix_id,...
    %     'sig_id',sig_id,'omg_state_id',omg_state_id,'omg_matrix_id',omg_matrix_id);
    
    [initvals,params.smoothed_vals,stoch_vol_small]=load_initial_conditions();
    % add priors for initial conditions
    return
end

% params: vector of parameters
% y : data
% X : matrix of regressors

% hard code the markov state to 1
%--------------------------------
% this may change later if we start markov-switching the parameters
istate=1;

%****************************
% PART I: DRAW THE VARIABLES
%****************************

% round 1: draw (A_t|)
%---------------------------------
if obj.time_varying_parameters(1)
    THETA_A_param=params(pp.theta_a_state_id,istate);
    if obj.random_walk_parameters(1)
        na=numel(THETA_A_param);
        A_param=zeros(na,1);
        RHO_A_param=ones(na,1);
    else
        A_param=params(pp.a_state_id,istate);
        RHO_A_param=params(pp.rho_a_state_id,istate);
    end
    % rewrite the system is state-space form with mean adjustments
    % we follow the notation in Durbin-Koopman
    %-------------------------------------------------------------
    c=(1-RHO_A_param).*A_param; % constant in the state equation
    d=[]; % constant in the measurement equation
    T={sparse(diag(RHO_A_param))}; % state matrix (autoregressive)
    R={sparse(diag(THETA_A_param))}; % state matrix (shocks)
    r=size(R{1},2);
    Q={eye(r)}; % covariance of shocks in state equation
    if obj.time_varying_parameters(3)
        OMGt=smoothed_vals(endo_ids.omg_state_id,:);
    else
        OMGt=params(pp.omg_state_id,istate);
    end
    if obj.time_varying_parameters(2)
        SIGt=smoothed_vals(endo_ids.sig_id,:);
    else
        SIGt=sqrt(params(pp.sig_state_id,istate));
    end
    % build Ht=(OMGt*SIGt)*(OMGt*SIGt)'
    H=build_H(OMGt,SIGt);% covariance of measurement errors
    % initial conditions for simulation smoothing
    %--------------------------------------------
    a1=initvals.At.a1; % a1=zeros(m,1);
    P1=initvals.At.P1; % P1=10*eye(m);
    % generate a simulation for At
    %-----------------------------
    sim_smooth=simulation_smoother(y,Z,T,R,H,Q,P1,a1,c,d);
    At=sim_smooth.alpha;
    % store immediately
    %------------------
    smoothed_vals(endo_ids.a_state_id,:)=At;
    %     At = stochvol_tools.draw_B(y,X,OMG,SIGt,diag(sig_b.^2),priors);
    % At includes observations from 0 to T : the inital conditions have to
    % be stored!!!
    % check that the block exogeneity is respected otherwise do something with
    % the carter_kohn algorithm!!!!!
else
    % estimate parameters A and store immediately
    %--------------------------------------------
    % estimate A with a GLS
    A=estimate_A();
    params(pp.a_state_id,istate)=A;
    At=A;
end

% round 2: draw (OMGt|)
%-------------------------------
% SIG_deflate_forecast_error
fe=sigma_deflate_forecast_error();
if obj.time_varying_parameters(3)
    THETA_OMG_param=params(pp.theta_omg_state_id,istate);
    if obj.random_walk_parameters(3)
        na=numel(THETA_OMG_param);
        OMG_param=zeros(na,1);
        RHO_OMG_param=ones(na,1);
    else
        OMG_param=params(pp.omg_state_id,istate);
        RHO_OMG_param=params(pp.rho_omg_state_id,istate);
    end
    % initial conditions for simulation smoothing
    %--------------------------------------------
    a1=initvals.OMGt.a1;
    P1=initvals.OMGt.P1;
    OMGt=estimate_OMG(fe,OMG_param,RHO_OMG_param,THETA_OMG_param,a1,P1);
    % store immediately
    %------------------
    smoothed_vals(endo_ids.omg_state_id,:)=OMGt;
else
    % estimate parameters OMG and store immediately
    %----------------------------------------------
    OMG=estimate_OMG(fe);
    params(pp.omg_state_id,istate)=OMG;
    OMGt=OMG;
end

% round 3: draw (SIGt|)
%-------------------------------
if obj.time_varying_parameters(2)
    if obj.random_walk_parameters(2)
        SIG2_param=ones(n,1);
        RHO_SIG_param=ones(n,1);
    else
        SIG2_param=params(pp.sig_state_id,istate);
        RHO_SIG_param=params(pp.rho_sig_state_id,istate);
    end
    THETA_SIG_param=params(pp.theta_sig_state_id,istate);
    % OMG_deflate_forecast_error
    fe=omega_deflate_forecast_error();
    lfe2=log(fe.^2+stoch_vol_small);
    % system is now: log(fe_t^2)=1*log(sig_t^2)+0.5*log(ey_t^2)
    %                log(sig_t^2)=(1-rho_sig)*log(sig^2)+rho_sig*log(sig_{t-1}^2)+theta_sig*upsil_sig_t
    % rewrite the system is state-space form with mean adjustments
    % we follow the notation in Durbin-Koopman
    %-------------------------------------------------------------
    c=(1-RHO_SIG_param).*log(SIG2_param); % constant in the state equation
    d=[]; % constant in the measurement equation
    T={sparse(diag(RHO_SIG_param))}; % state matrix (autoregressive)
    R={sparse(diag(THETA_SIG_param))}; % state matrix (shocks)
    Q={eye(n)}; % covariance of shocks in state equation
    ZZ={eye(n)};
    H={0.5^2*eye(n)};% covariance of measurement errors
    % initial conditions for simulation smoothing
    %--------------------------------------------
    a1=initvals.SIGt.a1; % a1=zeros(n,1);
    P1=initvals.SIGt.P1; % P1=10*eye(n);
    % generate a simulation
    %----------------------
    sim_smooth=simulation_smoother(lfe2,ZZ,T,R,H,Q,P1,a1,c,d);
    SIGt=sqrt(exp(sim_smooth.alpha));
    % store immediately
    %------------------
    smoothed_vals(endo_ids.sig_id,:)=SIGt;
    if any(isnan(SIGt(:)))
        keyboard
    end
else
    % estimate parameters OMG and store immediately
    %----------------------------------------------
    % estimate SIG with an OLS with a prior
    SIG2=estimate_SIG2();
    params(pp.sig_state_id,istate)=SIG2;
    SIGt=sqrt(SIG2);
end

%*****************************
% PART II: DRAW THE PARAMETERS
%*****************************
% round 4: draw (A,RHO_A,THETA_A|)
%---------------------------------
if obj.time_varying_parameters(1)
    [A_param,RHO_A_param,THETA_A_param]=single_equation_ols(At,obj.random_walk_parameters(1));
    if ~obj.random_walk_parameters(1)
        params(pp.a_state_id,istate)=A_param;
        params(pp.rho_a_state_id,istate)=RHO_A_param;
    end
    params(pp.theta_a_state_id,istate)=THETA_A_param;
end

% round 5. draw (SIG2,RHO_SIG and THETA_SIG|)
%-------------------------------------------
if obj.time_varying_parameters(2)
    [SIG2_param,RHO_SIG_param,THETA_SIG_param]=single_equation_ols(log(SIGt.^2),obj.random_walk_parameters(2));
    if ~obj.random_walk_parameters(2)
        params(pp.sig_state_id,istate)=SIG2_param;
        params(pp.rho_sig_state_id,istate)=RHO_SIG_param;
    end
    params(pp.theta_sig_state_id,istate)=THETA_SIG_param;
end

% round 6: draw (OMG,RHO_OMG,THETA_OMG|)
%---------------------------------------
if obj.time_varying_parameters(3)
    [OMG_param,RHO_OMG_param,THETA_OMG_param] = single_equation_ols(OMGt,obj.random_walk_parameters(3));
    if ~obj.random_walk_parameters(3)
        params(pp.omg_state_id,istate)=OMG_param;
        params(pp.rho_omg_state_id,istate)=RHO_OMG_param;
    end
    params(pp.theta_omg_state_id,istate)=THETA_OMG_param;
end

    function [initvals,smoothed_vals,stoch_vol_small]=load_initial_conditions()
        % Singularity adjustment factor (offset constant)
        %------------------------------------------------
        % the 'offset constant' prevents the case where fe is zero (in this
        % case, lfe2 is -Infinity). Primicieri(2005) hard-coded it to be
        % 0.001 
        bols=y/X;
        resid=y-bols*X;
        stoch_vol_small=min(sum(resid.^2,2)/smpl)/10000;

        % initial values for the smoothing simulations
        %---------------------------------------------
        smoothed_vals=ones(obj.endogenous.number(end),smpl);
        
        % initial conditions for the kalman filter
        %-----------------------------------------
        initvals=struct();
        if obj.time_varying_parameters(1)
            na_=numel(endo_ids.a_state_id);
            initvals.At.a1=zeros(na_,1);
            if obj.random_walk_parameters(1)
            else
            end
            initvals.At.P1=10*eye(na_);
        end
        if obj.time_varying_parameters(2)
            nsig_=numel(endo_ids.sig_id);
            initvals.SIGt.a1=zeros(nsig_,1);
            if obj.random_walk_parameters(2)
            else
            end
            initvals.SIGt.P1=10*eye(nsig_);
        end
        if obj.time_varying_parameters(3)
            nomg_=numel(endo_ids.omg_state_id);
            initvals.OMGt.a1=zeros(nomg_,1);
            if obj.random_walk_parameters(3)
            else
            end
            initvals.OMGt.P1=10*eye(nomg_);
        end
    end

    function [const,arcoefs,stdev]=single_equation_ols(data,random_walk)
        [nv,bigt]=size(data);
        const=zeros(nv,1);
        arcoefs=ones(nv,1);
        if random_walk
            resid=data(:,2:end)-data(:,1:end-1);
            stdev=sqrt(sum(resid.^2,2)/(bigt-1));
        else
            stdev=zeros(nv,1);
            constant_flag=true;
            nlags=1;
            xdata=[];
            for v=1:nv
                [results,stats]=vartools.ols(data(v,:),xdata,nlags,constant_flag);
                arcoefs(v)=results.A(1);
                const(v)=results.A0;
                SIGols=results.SIGols;
                stdev(v)=sqrt(SIGols);
            end
            const=const./(1-arcoefs);
        end
    end

    function Z=build_Z()
        Z=cell(1,smpl);
        In=eye(n);
        vec_order=cell2mat(obj.parameters_positions.a{2}');
        vec_order=vec_order(:,3);
        for t=1:smpl
            Zt=kron(X(:,t)',In);
            % re-order the columns and pick only the ones corresponding to
            % the non-zero elements in the A+ matrix
            %-------------------------------------------------------------
            Z{t}=sparse(Zt(:,vec_order));
        end
    end

    function H=build_H(OMG,SIG)
        smax=size(SIG,2);
        omax=size(OMG,2);
        bigt=max(smax,omax);
        H=cell(1,bigt);
        for t=1:bigt
            if smax>=t
                SIGt_=diag(SIG(:,t));
            end
            if omax>=t
                OMGt_=build_OMG(OMG(:,t));
            end
            H{t}=OMGt_*SIGt_;
            H{t}=H{t}*H{t}';
        end
    end

    function fe=sigma_deflate_forecast_error()
        fe=forecast_error();
        fe=bsxfun(@rdivide,fe,SIGt);
    end

    function fe=omega_deflate_forecast_error()
        fe=forecast_error();
        if size(OMGt,2)>1 % ||obj.time_varying_parameters(3)
            for t=1:smpl
                fe(:,t)=build_OMG(OMGt(:,t))\fe(:,t);
            end
        else
            fe=build_OMG(OMGt)\fe;
        end
    end

    function fe=forecast_error()
        fe=y;
        if size(At,2)>1 % ||obj.time_varying_parameters(1)
            for t=1:smpl
                fe(:,t)=y(:,t)-build_A(At(:,t))*X(:,t);
            end
        else
            fe=y-build_A(At)*X;
        end
    end

    function OMGt=estimate_OMG(fe,OMG_param,RHO_OMG_param,THETA_OMG_param,a1,P1)
        y1=fe(1,:);
        xdata=y1(ones(1,n),:);
        OMG_=eye(n);
        is_time_varying=nargin>1;
        if is_time_varying
            OMG_=repmat(OMG_,[1,1,smpl]);
            % rewrite the system is state-space form with mean adjustments
            % we follow the notation in Durbin-Koopman
            %-------------------------------------------------------------
            cc=build_OMG((1-RHO_OMG_param).*OMG_param); % constant in the state equation
            di=[]; % constant in the measurement equation
            TT=build_OMG(RHO_OMG_param); % state matrix (autoregressive)
            RR=build_OMG(THETA_OMG_param); % state matrix (shocks)
            Hi={1};% covariance of measurement errors
            a1=build_OMG(a1);
            P1=build_OMG(diag(P1));
            dim1Dist=1;
            for irow=2:n
                a1i=a1(irow,1:irow-1); a1i=a1i(:);
                P1i=P1(irow,1:irow-1); P1i=diag(P1i(:));
                ci=cc(irow,1:irow-1); ci=ci(:);
                Ti={diag(TT(irow,1:irow-1))};
                Ri={diag(RR(irow,1:irow-1))};
                rr=size(Ri{1},2);
                Qi={eye(rr)}; % covariance of shocks in state equation
                tmp=xdata(1:irow-1,:);
                dim2Dist=(irow-1)*ones(1,smpl);
                Zi=mat2cell(tmp(:)',dim1Dist,dim2Dist);
                ytilde=fe(irow,:);
                % generate a simulation for At
                %-----------------------------
                sim_smooth=simulation_smoother(ytilde,Zi,Ti,Ri,Hi,Qi,P1i,a1i,ci,di);
                OMG_(irow,1:irow-1,:)=sim_smooth.alpha;
                for t=1:smpl
                    yhat=Zi{t}*sim_smooth.alpha(:,t);
                    xdata(irow,t)=fe(irow,t)-yhat;
                end
                xdata(irow,:)=ytilde-xdata(irow,:);
            end
            OMGt=vectorize_OMG(OMG_(:,:,1));
            OMGt(:,ones(1,smpl));
            for t=2:smpl
                OMGt(:,t)=vectorize_OMG(OMG_(:,:,t));
            end
        else
            constant_flag=false;
            nlags=0;
            for irow=2:n
                ytilde=fe(irow,:);
                [results,stats]=vartools.ols(ytilde,xdata(1:irow-1,:),nlags,constant_flag);
                if constant_flag
                    arcoefs_=results.A(1:end-1);
                else
                    arcoefs_=results.A(1:end);
                end
                res=results.residuals;
                
                OMG_(irow,1:irow-1)=arcoefs_(:)';
                xdata(irow,:)=res;
            end
            OMGt=vectorize_OMG(OMG_);
        end
    end
end
