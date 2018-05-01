function [loglik,Incr,retcode,Filters]=switching_divided_difference_filter(syst,y,U,z,options)
% switching_divided_difference_filter - filter for nonlinear regime-swiching models
%
% ::
%
%
%   [loglik,Incr,retcode,Filters]=switching_divided_difference_filter(...
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
%    - **y** [matrix]: ny x T matrix of data
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
% Example:
%
%    See also:

% hard-coded for now
type='stirling';

pi0=syst.PAI00;

x0=syst.a;

h = numel(x0);

Qfunc=syst.Qfunc;

ff=syst.ff;

kf_householder_chol=options.kf_householder_chol;

S0=cell_householder(syst.P);

Sw=cell_householder(syst.H);

Sv=cell_householder(syst.SIGeta);

store_filters=options.kf_filtering_level;
% smoothing not ready
%---------------------
store_filters=min(store_filters,2);

if store_filters
    
    nsteps=options.kf_nsteps;
    
else
    % do not do multi-step forecasting during estimation
    nsteps=1;
    
end

kalman_tol=options.kf_tol;

% first iterate to consider in the computation of the likelihood.
%-------------------------------------------------------------------
first=syst.start;

nx = size(S0{1},1);

nv = size(Sv{1},2); % number of structural shocks

% zeros shocks
%----------------
horizon=syst.horizon;

shocks0=zeros(nv,horizon);

% deterministic variables only if rows present
%----------------------------------------------
tmax_u=size(U,2)*(size(U,1)>0);

if tmax_u
    % update the state (x_{1|0} with the deterministic variables
    %------------------------------------------------------------
    Ut=U(:,0+1);
    
    for splus=1:h
        
        [x0{splus}]=ff(splus,x0{splus},shocks0,Ut);
    
    end
    
end

[ny0,T] = size(y);

delta = sqrt(3);

scal0 = (delta^2-nx-nv)/(delta^2);

scal1 = 1/(2*delta^2);

scal2 = 1/(2*delta);

scal3 = sqrt(delta^2-1)*scal1;

% initialize various fields
%---------------------------
Sbarx = S0;

xbar_tt1 = x0;

x_hat = cell(h,1);

Shat_x = cell(h,1);

% for the forecasting step
Sxx_1 = zeros(nx);

Sxx_2 = zeros(nx);

Sxv_1 = zeros(nx,nv);

Sxv_2 = zeros(nx,nv);

Incr = zeros(T,1);

loglik=[];

att_all=x0{1}*pi0(1);

for rt=2:h
    
    att_all=att_all+x0{rt}*pi0(rt);
    
end

Q=Qfunc(att_all);

pi0=transpose(Q)*pi0;

[Filters,K_store,iF_store,v_store]=utils.filtering.initialize_storage(x0,S0,...
    pi0,Q,ny0,nv,horizon,nsteps,T,store_filters);

pi_tt1 = pi0;

retcode=0;

for t = 1:T
    % data and indices for observed variables at time t
    %--------------------------------------------------
    [ny,occur,obsOccur]=z(t);
    
    yt=y(occur,t);
    
    log_f01 = nan(h,1);
    
    for rt = 1:h
        
        % forecast of observables
        %-------------------------
        y_tt1=xbar_tt1{rt}(obsOccur); %<-- y_tt1 = z*xbar_tt1{rt};
        
        % forecast error
        %-----------------
        v = yt - y_tt1;
        
        % covariance of forecast error
        %------------------------------
        Syxbar_1 = Sbarx{rt}(obsOccur,:); %<----Syxbar_1 = z*Sbarx{rt};
        
        Syw_1 = Sw{rt};
        
        Pxy = Sbarx{rt}*Syxbar_1';
        
        Sy = utils.cov.householder([Syxbar_1,Syw_1],kf_householder_chol);
        
        F = Sy*Sy';
        
        detF=det(F);
        
        failed=detF<=0;
        
        if ~failed
            % inverse of covariance of forecast error
            %----------------------------------------
            Finv = eye(ny)/F;
            
            failed=any(isnan(Finv(:)));
            
        end
        
        if failed
            
            retcode=305;
            
            return
            
        end
        
        % Kalman gain
        %-------------
        K = Pxy*Finv;
        
        % state update
        %--------------
        x_hat{rt} = xbar_tt1{rt} + K*v;
        
        % covariance of state update
        %----------------------------
        if isempty(Syw_1)
            
            KSyw_1=[];
            
        else
            
            KSyw_1=K*Syw_1;
            
        end
        
        Shat_x{rt} = utils.cov.householder([Sbarx{rt} - K*Syxbar_1, KSyw_1],...
            kf_householder_chol);
        
        % contribution to the likelihood
        %--------------------------------
        log_f01(rt)=-0.5*(...
            log((2*pi)^ny*detF)+...
            v'*Finv*v...
            );
        
    end
    
    [Incr(t),PAI01_tt,retcode]=switch_like_exp_facility(pi_tt1,log_f01,kalman_tol);
    
    if retcode
        
        return
        
    end
    
    % update of probabilities
    %-------------------------
    pi_tt = sum(PAI01_tt,2);
    
    if store_filters>1
        
        store_updates();
        
    end

    % update of transition matrix
    %-----------------------------
    if h>1
        
        xhat_all=x_hat{1}*pi_tt(1);
        
        for rt=2:h
            
            xhat_all=xhat_all+x_hat{rt}*pi_tt(rt);
            
        end
        
        [Qtt,retcode]=Qfunc(xhat_all);
        
        if retcode
            
            return
            
        end
        
        % forecast of probabilities
        %---------------------------
        pi_tt1 = Qtt'*pi_tt;
        
    else
        
        Qtt=1;
        
        pi_tt1 =1;
        
    end
    
    for rt1 = 1:h
        
        % collapsing
        %------------
        xc_t1 = 0;
        
        Sxc_t1 = 0;
        
        Scv_t1 = 0;
        
        for rt = 1:h
            
            if h==1
                
                pi_rtrt1=1;
                
            else
                
                pi_rtrt1 = Qtt(rt,rt1)*pi_tt(rt)/pi_tt1(rt1);
                
            end
            
            xc_t1 = xc_t1 + pi_rtrt1*x_hat{rt};
            
            Sxc_t1 = Sxc_t1 + pi_rtrt1*Shat_x{rt};
            
            Scv_t1 = Scv_t1 + pi_rtrt1*Sv{rt};
            
        end
        
        if t+1<=tmax_u
            
            Ut=U(:,t+1);
            
        else
            
            Ut=[];
            
        end
        
        [Sbarx{rt1}, xbar_tt1{rt1}] = forecast_cov(xc_t1,Sxc_t1,Scv_t1,rt1);
                
    end
    
    if store_filters>0
        
        store_predictions()
        
        if store_filters>2
            
%             R_store(t).R=Rt;

        end
        
    end

   if options.debug
       
       disp(t)
       
   end
   
end

loglik = sum(Incr(first:end));

    function [Sbarx, xbar] = forecast_cov(x,Sx,Sv,rt1)
        
        xbar = 0;
        
        f0=ff(rt1,x,shocks0,Ut);
        
        for p = 1:nx
            
            fplus=ff(rt1,x + delta*Sx(:,p),shocks0,Ut);
            
            fminus=ff(rt1,x - delta*Sx(:,p),shocks0,Ut);
            
            Sxx_1(:,p) = ( fplus - fminus );
            
            Sxx_2(:,p) = ( fplus + fminus - 2*f0);
            
            xbar = xbar + fplus + fminus;
            
        end
        
        Sxx_1=scal2*Sxx_1;
        
        Sxx_2=scal3*Sxx_2;
        
        for p = 1:nv
            
            fplus=ff(rt1,x,delta*Sv(:,p),Ut);
            
            fminus=ff(rt1,x, - delta*Sv(:,p),Ut);
            
            Sxv_1(:,p) = ( fplus - fminus );
            
            Sxv_2(:,p) = ( fplus + fminus - 2*f0 );
            
            xbar = xbar + fplus + fminus;
            
        end
        
        Sxv_1=scal2*Sxv_1;
        
        Sxv_2=scal3*Sxv_2;
        
        Sbarx = utils.cov.householder([Sxx_1,Sxv_1,Sxx_2,Sxv_2],kf_householder_chol);
        
        if strcmp(type,'stirling')
            
            xbar = scal0*f0 +  scal1*xbar;
            
        elseif strcmp(type,'none')
            
            xbar = f0;
            
        end
        
    end

    function store_predictions()
        
        Filters.PAI(:,t+1)=pi_tt1;
        
        Filters.Q(:,:,t+1)=Qtt;
        
        for splus_=1:h
            
            Filters.a{splus_}(:,1,t+1)=xbar_tt1{splus_};
            
            Filters.P{splus_}(:,:,t+1)=Sbarx{splus_}*Sbarx{splus_}.';
            
            for istep_=2:nsteps
                % this assumes that we stay in the same state and we know
                % we will stay. The more general case where we can jump to
                % another state is left to the forecasting routine. Here we
                % just do the conditional... not the full expectation...
                if t+istep_-1<=tmax_u
                    
                    Utplus=U(:,t+istep_-1);
                    
                else
                    
                    Utplus=[];
                    
                end
                
                [Filters.a{splus_}(:,istep_,t+1),~,rcode]=ff(splus_,...
                    Filters.a{splus_}(:,istep_-1,t+1),shocks0,Utplus);

                if rcode
                    % do not exit completely...
                    break
                    
                end
                
            end
            
        end
        
    end

    function store_updates()
        
        Filters.PAItt(:,t)=pi_tt;
        
        for st_=1:h
            
            Filters.att{st_}(:,1,t)=x_hat{st_};
            
            Filters.Ptt{st_}(:,:,t)=Shat_x{st_}*Shat_x{st_}.';
            
            if store_filters>2
                
                K_store{st_}(:,occur,t)=K(:,occur,st_);
                
                iF_store{st_}(occur,occur,t)=iF{st_};
                
                v_store{st_}(occur,t)=v{st_};
                
            end
            
        end
        
    end

    function S=cell_householder(P)
        
        S=cell(1,h);
        
        for reg=1:h
            
            if ~isempty(P{reg})
                
                [UU,SS,VV] = svd(P{reg});
                
                C=utils.cov.householder(UU*diag(sqrt(diag(SS)))*VV.',...
                    kf_householder_chol);
                
                S{reg}=C; % <-- S{rt}=chol(P{rt},'lower');
                
            end
            
        end
        
    end

end

