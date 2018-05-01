function f=missing_observations_kalman_filter(y,T,R,Z,H,Q,init,dy,ca,level,nstep)
% MISSING_OBSERVATIONS_KALMAN_FILTER -- Kalman filter for
% non-regime-switching models
%
% ::
%
%
%   f=MISSING_OBSERVATIONS_KALMAN_FILTER(y,T,R,Z,H,Q,init,dy,ca,level)
%
% Args:
%
%      - **y** [p x n matrix]: data
%
%      - **T** [m x m x nT]: state matrix, see below
%
%      - **R** [m x r x nR]: state matrix, see below
%
%      - **Z** [p x m||logical]: state matrix or selector, see below
%
%      - **H** [p x p x nH]: Covariances of measurement errors
%
%      - **Q** [r x r x nQ]: Covariances of shocks
%
%      - **init** [struct]:
%
%      - **dy** [p x ndy]:
%
%      - **ca** [m x nca]:
%
%      - **level** [0|1|2|{3}]: controls the level output produced.
%          if level==0, only the likelihood is returned
%          if level==1, filters/forecasts are returned in addition to level 0
%          if level==2, updates are returned in addition to level 1
%          if level==3, smoothed values are returned in addition to level 2
%
%      - **nstep** [integer|{1}]: Number of forecast steps
%
% Returns:
%    :
%
%      - **f** [struct]: output structure
%          - **retcode** [0|305]: 305 if computation of the likelihood breaks
%          downs
%          - **incr** [vector]: likelihood in each period
%          - **log_lik** [scalar]: log likelihood
%          - **a** [m x n+1 x nstep]: filtered stated (available if level >0)
%          - **P** [m x m x n+1 x nstep]: Covariance matrix of filtered state
%          (available if level >0)
%          - **att** [m x n]: Updated state vector  (available if level >1)
%          - **Ptt** [m x m x n]: Covariance matrix for updates (available if
%          level >1)
%          - **iF** [p x p x n]: Inverses of covariance matrix of forecast
%          errors
%          - **v** [p x n]: Forecast errors
%          - **K** [m x p x n]: Kalman gains
%          - **alpha** [m x n]: smoothed state vector
%          - **V** [m x m x n]: Covariance matrix of smoothed state vector
%          - **r** [m x n+1]: quantity useful for the efficient computation of
%          the state vector
%          - **eta** [r x n]: smoothed shocks
%          - **epsilon** [p x n]: smoothed measurement errors
%
% Note:
%
%    The system is of the form
%      alpha(t+1) = T*alpha(t) + ca(t) + R(t)*eta(t),  eta ~ N(0, Q)
%      y(t) =   Z*alpha(t) + dy(t) + epsilon(t),  epsilon ~ N(0,H)
%
% Example:
%
%    See also:

if nargin < 11
    
    nstep=[];
    
    if nargin<10
        
        level=[];
        
        if nargin<9
            
            ca=[];
            
            if nargin<8
                
                dy=[];
                
                if nargin<7
                    
                    init=[];
                    
                end
                
            end
            
        end
        
    end
    
end

if isempty(nstep), nstep=1; end

if isempty(level),level=3; end

if ~any(level==[0,1,2,3])
    
    error('level must be 0, 1, 2, or 3')
    
end

ndy=size(dy,2); nca=size(ca,2);

[p,n]=size(y);

good=~isnan(y);

m=size(T,1);

r=size(R,2);

[sp,Z,Z_is_selector]=obsw.time_series_length(Z,T,H,Q,R);

if Z_is_selector
    
    Im=eye(m);
    
    Zmat=Im(Z,:);
    
end

[a1,P1]=kalman_initialization(init);

f=output_initialization(a1,P1);

the_filter(a1,P1)

if level>2 && ~f.retcode
    
    the_smoother()
    
end

    function the_smoother()
        
        rt=zeros(m,1);
        
        Nt=zeros(m);
        
        for t=n:-1:1
            
            if Z_is_selector
                
                Zt=Zmat;
                
            else
                
                Zt=Z(:,:,sp.Z(t));
                
            end
            
            % structural shocks [begining or end of period?]
            f.eta(:,t)=Q(:,:,sp.Q(t))*R(:,:,sp.R(t)).'*rt;
            
            Tt=T(:,:,sp.T(t));
            
            Pt=f.P(:,:,t,1);
            
            Kt=f.K(:,:,t); % Tt*Pt*Zt.'*f.iF(:,:,t);
            
            Lt=Tt-Kt*Zt;
            
            ut=f.iF(:,:,t)*f.v(:,t)-Kt.'*rt;
            
            rt=Zt.'*f.iF(:,:,t)*f.v(:,t)+Lt.'*rt;
            % = Zt'*ut+Tt'*rt
            % = Zt'*(f.iF(:,:,t)*f.v(:,t)-Kt.'*rt)+Tt'*rt
            
            f.r(:,t)=rt;
            
            Nt=Zt.'*f.iF(:,:,t)*Zt+Lt.'*Nt*Lt;
            
            f.alpha(:,t)=f.a(:,t,1)+Pt*rt;
            
            f.V(:,:,t)=Pt-Pt*Nt*Pt;
            
            % measurement errors
            f.epsilon(:,t)=H(:,:,sp.H(t))*ut;
            
            % %             % structural shocks [begining or end of period?]
            % %             f.eta(:,t)=Q(:,:,sp.Q(t))*R(:,:,sp.R(t)).'*rt;
            
        end
        
    end

    function the_filter(a,P)
        
        l2pi=log(2*pi);
        
        for t=1:n
            
            Tt=T(:,:,sp.T(t));
            
            okt=good(:,t);
            
            if Z_is_selector
                
                Zt=Z(okt);
                
                v=y(okt,t)-a(Zt)-dy(okt,tdy(t));
                
                F=P(Zt,Zt)+H(okt,okt,sp.H(t));
                
            else
                
                Zt=Z(okt,:);
                
                v=y(okt,t)-Zt*a-dy(okt,tdy(t));
                
                F=Zt*P*Zt.'+H(okt,okt,sp.H(t));
                
            end
            
            pstar=sum(okt);
            
            detF=det(F);
            
            failed=detF<=0;
            
            if ~failed
                
                iF=F\eye(pstar);
                
                failed=any(isnan(iF(:)));
                
            end
            
            if failed
                
                f.retcode=305;
                
                return
                
            end
            
            f.incr(t)=-0.5*(pstar*l2pi+log(detF)+v.'*iF*v);
            %             f.incr(t)=-pstar/2*log(2*pi)-1/2*(log(det(F))+v.'*iF*v);
            
            update()
            
            forecast()
            
        end
        
        f.log_lik=sum(f.incr);
        
        function update()
            
            if Z_is_selector
                
                PZiF=P(:,Zt)*iF;
                
                Kt=Tt*PZiF; % Ref. page 85
                
                a=a+PZiF*v;
                
                P=P-PZiF*P(Zt,:);
                
            else
                
                PZiF=P*Zt.'*iF;
                
                Kt=Tt*PZiF; % Ref. page 85
                
                a=a+PZiF*v;
                
                P=P-PZiF*Zt*P;
                
            end
            
            if level>1
                
                f.att(:,t)=a;
                
                f.Ptt(:,:,t)=P;
                
                if level>2
                    
                    f.K(:,okt,t)=Kt;
                    
                    f.iF(okt,okt,t)=iF;
                    
                    f.v(okt,t)=v;
                    
                end
                
            end
            
        end
        
        function forecast()
            
            Rt=R(:,:,sp.R(t));
            
            a=Tt*a+ca(:,tca(t));
            
            P=Tt*P*Tt.'+Rt*Q(:,:,sp.Q(t))*Rt.';
            
            if level>0
                
                f.a(:,t+1,1)=a;
                
                f.P(:,:,t+1,1)=P;
                
                for ii=2:nstep
                    
                    [f.a(:,t+1,ii),f.P(:,:,t+1,ii)]=...
                        hairy(f.a(:,t+1,ii-1),f.P(:,:,t+1,ii-1));
                    
                end
                
            end
            
            function [a,P]=hairy(a,P)
                
                tau=t+ii-1;
                
                Tt2=T(:,:,sp.T(tau));
                
                Rt2=R(:,:,sp.R(tau));
                
                a=Tt2*a+ca(:,tca(tau));
                
                P=Tt2*P*Tt2.'+Rt2*Q(:,:,sp.Q(tau))*Rt2.';
                
            end
            
        end
        
    end

    function a=tdy(t),a=min(t,ndy); end

    function a=tca(t),a=min(t,nca); end

    function [a1,P1]=kalman_initialization(init)
        
        if isempty(ca),ca=zeros(m,1); nca=1; end
        
        if isempty(dy),dy=zeros(p,1); ndy=1;end
        
        if isempty(H),H=zeros(p); end
        
        if isempty(init)
            
            RQR=R(:,:,1)*Q(:,:,1)*R(:,:,1).';
            
            [P1,retcode]=doubling_solve(T(:,:,1),T(:,:,1).',RQR);
            
            if retcode
                
                error(decipher(retcode))
                
            end
            
            init=struct('a',zeros(m,1),'P',P1);
            
        end
        
        a1=init.a;
        
        P1=init.P;
        
    end

    function f=output_initialization(a1,P1)
        
        f=struct();
        
        f.incr=nan(n,1);
        
        f.log_lik=nan;
        
        f.retcode=0;
        
        if level>0
            
            f.a=zeros(m,n+1,nstep); f.a(:,1,1)=a1;
            
            f.P=zeros(m,m,n+1,nstep); f.P(:,:,1,1)=P1;
            
            for istep=2:nstep
                
                f.a(:,1,istep)=a1;
                
                f.P(:,:,1,istep)=P1;
                
            end
            
            if level>1
                
                f.att=zeros(m,n);
                
                f.Ptt=zeros(m,m,n);
                
                if level>2
                    
                    f.iF=zeros(p,p,n);
                    
                    f.v=zeros(p,n);
                    
                    f.K=zeros(m,p,n);
                    
                    f.alpha=zeros(m,n);
                    
                    f.V=zeros(m,m,n);
                    
                    f.r=zeros(m,n+1);
                    
                    f.eta=zeros(r,n);
                    
                    f.epsilon=zeros(p,n);
                    
                end
                
            end
            
        end
        
    end

end