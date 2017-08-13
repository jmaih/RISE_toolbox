function f=missing_observations_kalman_filter(y,T,R,Z,H,Q,init,dy,ca)

%   alpha(t+1) = T*alpha(t) + ca(t) + R(t)*eta(t),  eta ~ N(0, Q)
%   y(t) =   Z*alpha(t) + epsilon(t),  epsilon ~ N(0,H)
%
if nargin<9
    
    ca=[];
    
    if nargin<8
        
        dy=[];
        
        if nargin<7
            
            init=[];
            
        end
        
    end
    
end

ndy=size(dy,2); nca=size(ca,2); 

[p,n]=size(y);

good=~isnan(y);

m=size(T,1);

[sp,Z,Z_is_selector]=obsw.time_series_length(Z,T,H,Q,R);

if Z_is_selector
    
    Im=eye(m);
    
    Zmat=Im(Z,:);
    
end

[a1,P1]=kalman_initialization(init);

f=output_initialization(a1,P1);

the_filter(f.a(:,1),f.P(:,:,1))

the_smoother()

    function the_smoother()
                
        rt=zeros(m,1);
        
        Nt=zeros(m);
        
        f.alpha=zeros(m,n);
        
        f.V=zeros(m,m,n);
        
%         f.N=zeros(m,m,n+1);
        
%         f.C=zeros(m,p,n);
        
%         f.r=zeros(m,n+1);
        
        for t=n:-1:1
            
            if Z_is_selector
                
                Zt=Zmat;
                
            else
                
                Zt=Z(:,:,sp.Z(t));
                
            end
            
            % structural shocks [begining or end of period?]
            f.eta(:,t)=Q(:,:,sp.Q(t))*R(:,:,sp.R(t)).'*rt;
            
            Tt=T(:,:,sp.T(t));
            
            Pt=f.P(:,:,t);
            
            Kt=f.K(:,:,t); % Tt*Pt*Zt.'*f.iF(:,:,t);
            
%             Dt=f.iF(:,:,t)+Kt.'*Nt*Kt;
            
%             f.C(:,:,t)=Zt.'*Dt-Tt.'*Nt*Kt;
            
            Lt=Tt-Kt*Zt;
            
            ut=f.iF(:,:,t)*f.v(:,t)-Kt.'*rt;
            
            rt=Zt.'*f.iF(:,:,t)*f.v(:,t)+Lt.'*rt; 
            % = Zt'*ut+Tt'*rt
            % = Zt'*(f.iF(:,:,t)*f.v(:,t)-Kt.'*rt)+Tt'*rt
            
%             f.r(:,t)=rt;
            
            Nt=Zt.'*f.iF(:,:,t)*Zt+Lt.'*Nt*Lt;
            
%             f.N(:,:,t)=Nt;
            
            f.alpha(:,t)=f.a(:,t)+Pt*rt;
            
            f.V(:,:,t)=Pt-Pt*Nt*Pt;
            
            % measurement errors
            f.epsilon(:,t)=H(:,:,sp.H(t))*ut;
            
% %             % structural shocks [begining or end of period?]
% %             f.eta(:,t)=Q(:,:,sp.Q(t))*R(:,:,sp.R(t)).'*rt;
        
        end
        
    end

    function the_filter(a,P)
        
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
            
            iF=F\eye(pstar);
            
            update()
            
            forecast()
            
        end
        
        function update()
            
            if Z_is_selector
                
                f.K(:,okt,t)=Tt*P(:,Zt)*iF; % Ref. page 85
                
                a=a+P(:,Zt)*iF*v;
                
                P=P-P(:,Zt)*iF*P(Zt,:);
                
            else
                
                f.K(:,okt,t)=Tt*P*Zt.'*iF; % Ref. page 85
                
                a=a+P*Zt.'*iF*v;
                
                P=P-P*Zt.'*iF*Zt*P;
                
            end
            
            f.att(:,t)=a;
            
            f.Ptt(:,:,t)=P;
            
            f.iF(okt,okt,t)=iF;
            
            f.v(okt,t)=v;
            
        end
        
        function forecast()
            
            Rt=R(:,:,sp.R(t));
            
            a=Tt*a+ca(:,tca(t));
            
            P=Tt*P*Tt.'+Rt*Q(:,:,sp.Q(t))*Rt.';
            
            f.a(:,t+1)=a;
            
            f.P(:,:,t+1)=P;
            
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
        
        f.a=zeros(m,n+1); f.a(:,1)=a1;
        
        f.P=zeros(m,m,n+1); f.P(:,:,1)=P1;
        
        f.att=zeros(m,n);
        
        f.Ptt=zeros(m,m,n);
        
        f.iF=zeros(p,p,n);
        
        f.v=zeros(p,n);
        
        f.K=zeros(m,p,n);
        
    end

end