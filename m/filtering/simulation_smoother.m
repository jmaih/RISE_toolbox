function [sim_smooth,smooth,filt]=simulation_smoother(y,Z,T,R,H,Q,P1,a1,c,d)
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

% the model 
% y_t=Z_t*alpha_t+d_t+e_t
% alpha_{t+1}=T_t*alpha_t+c_t+R_t*eta_t
%{
m=10;
p=3;
n=100;
r=5;
c=randn(m,n);
d=randn(p,n);
y=nan(p,n);
alpha=zeros(m,n+1);
Z=cell(1,n); T=cell(1,n); R=cell(1,n); H=cell(1,n); Q=cell(1,n);
for t=1:n
    Z{t}=randn(p,m);
    T{t}=randn(m);
    T{t}=T{t}/norm(T{t});
    R{t}=randn(m,r);
    H{t}=diag(rand(p,1));
    Q{t}=diag(rand(r,1));
    y(:,t)=Z{t}*alpha(:,t)+d(:,t)+sqrt(H{t})*randn(p,1);
    alpha(:,t+1)=T{t}*alpha(:,t)+c(:,t)+R{t}*sqrt(Q{t})*randn(r,1);
end
a1=zeros(m,1);
P1=lyapunov_equation(T{1},R{1}*Q{1}*R{1}');

[sim_smooth,smooth,filt]=simulation_smoother(y,Z,T,R,H,Q,P1,a1,c,d)

% Fitted observables 
%-------------------
yfitted_s=y;
for t=1:n
    yfitted_s(:,t)=Z{t}*smooth.alpha(:,t)+d(:,t);
end
% smoother
%---------
max(max(abs(smooth.e+yfitted_s-y)))
%}
if nargin<10
    
    d=[];
    
    if nargin<9
        
        c=[];
        
    end
    
end

[p,n]=size(y);

m=size(T{1},1);

r=size(R{1},2);

zmax=numel(Z);

tmax=numel(T);

rmax=numel(R);

qmax=numel(Q);

hmax=numel(H);

cmax=size(c,2);

dmax=size(d,2);

% filter
%--------
filt=do_filtering_step(y);

% smoothing step
%---------------
smooth=do_smoothing_step(filt);

sim_smooth=do_simulation_smoothing_step();

    function smooth=do_smoothing_step(filt)
        
        rt=zeros(m,1);
        
        alpha=zeros(m,n);
        
        e=zeros(p,n);
        
        eta=zeros(r,n);
        
        Im=eye(m);
        
        for t=n:-1:1
            % matrices
            %---------
            Zt=Z{min(t,zmax)};
            
            Tt=T{min(t,tmax)};
            
            Rt=R{min(t,rmax)};
            
            Qt=Q{min(t,qmax)};
            
            Ht=H{min(t,hmax)};
            
            Kt=Tt*filt.Kbar(:,:,t); % Lt=Tt-Tt*Kbar(:,:,t)*Zt;
            
            Lt=Tt*(Im-filt.Kbar(:,:,t)*Zt);
            
            ut=filt.iF(:,:,t)*filt.v(:,t)-Kt'*rt;
            
            e(:,t)=Ht*ut;
            
            eta(:,t)=Qt*Rt'*rt;
            
            rt=Zt'*filt.iF(:,:,t)*filt.v(:,t)+Lt'*rt;
            
            alpha(:,t)=filt.a(:,t)+filt.P(:,:,t)*rt;
            
        end
        
        smooth=struct('alpha',alpha,'e',e,'eta',eta);
        
    end

    function filt=do_filtering_step(y)
        % initialization
        %---------------
        P=zeros(m,m,n+1);
        
        P(:,:,1)=P1;
        
        a=zeros(m,n+1);
        
        a(:,1)=a1;
        
        v=zeros(p,n);
        
        iF=zeros(p,p,n);
        
        Kbar=zeros(m,p,n);
        
        Ip=eye(p);
        
        for t=1:n
            % matrices
            %---------
            if zmax>=t,Zt=Z{min(t,zmax)};end
            
            if tmax>=t,Tt=T{min(t,tmax)};end
            
            if rmax>=t,Rt=R{min(t,rmax)};end
            
            if qmax>=t,Qt=Q{min(t,qmax)};end
            
            if hmax>=t,Ht=H{min(t,hmax)};end
            
            % forecast error
            %---------------
            v(:,t)=y(:,t)-Zt*a(:,t);
            
            if dmax>=t
                
                v(:,t)=v(:,t)-d(:,t);
                
            end
            
            Ft=Zt*P(:,:,t)*Zt'+Ht;
            
            iF(:,:,t)=Ft\Ip;
            
            % updating step
            %---------------
            Kbar(:,:,t)=P(:,:,t)*Zt'*iF(:,:,t);
            
            att=a(:,t)+Kbar(:,:,t)*v(:,t);
            
            Ptt=P(:,:,t)-Kbar(:,:,t)*Zt*P(:,:,t);
            
            % filtering step
            %---------------
            a(:,t+1)=Tt*att;
            
            if cmax>=t
                
                a(:,t+1)=a(:,t+1)+c(:,t);
                
            end
            
            P(:,:,t+1)=Tt*Ptt*Tt'+Rt*Qt*Rt';
            
        end
        
        filt=struct('a',a,'P',P,'v',v,'iF',iF,'Kbar',Kbar);
        
    end

    function sim_smooth=do_simulation_smoothing_step()
        % draw shocks
        %------------
        CP1=efficient_chol(P1);
        
        alpha_plus=zeros(m,n+1);
        
        alpha_plus(:,1)=a1+CP1*randn(m,1);
        
        yplus=zeros(p,n);
        
        e_plus=zeros(p,n);
        
        eta_plus=zeros(r,n);
        
        for t=1:n
            % matrices
            %---------
            if zmax>=t,Zt=Z{min(t,zmax)};end
            
            if tmax>=t,Tt=T{min(t,tmax)};end
            
            if rmax>=t,Rt=R{min(t,rmax)};end
            
            if qmax>=t,Qt=Q{min(t,qmax)};CQt=efficient_chol(Qt);end
            
            if hmax>=t,Ht=H{min(t,hmax)};CHt=efficient_chol(Ht);end
            
            e_plus(:,t)=CHt*randn(p,1);
            
            eta_plus(:,t)=CQt*randn(r,1);
            
            yplus(:,t)=Zt*alpha_plus(:,t)+e_plus(:,t);
            
            alpha_plus(:,t+1)=Tt*alpha_plus(:,t)+Rt*eta_plus(:,t);
            
        end
        
        filt_=do_filtering_step(yplus);
        
        % smooth back
        %------------
        smooth_hat_plus=do_smoothing_step(filt_);
        
        % simulation smoothing variables
        %-------------------------------
        sim_smooth=struct();
        
        sim_smooth.alpha=alpha_plus(:,1:end-1)-smooth_hat_plus.alpha+smooth.alpha;
        
        sim_smooth.e=e_plus-smooth_hat_plus.e+smooth.e;
        
        sim_smooth.eta=eta_plus-smooth_hat_plus.eta+smooth.eta;
        
        function C=efficient_chol(X)
            
            dX=diag(X);
            
            if max(abs(diag(dX)-X))<1e-10
                
                C=diag(sqrt(dX));
                
            else
                
                C=chol(X,'lower');
                
            end
            
        end
        
    end

end
