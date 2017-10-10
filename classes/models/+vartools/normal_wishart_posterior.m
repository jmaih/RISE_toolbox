function [abar_out,SIGu,smpl]=normal_wishart_posterior(Z,SIGu,V,astar,y,linres)

if nargin < 6
    
    linres=[];
    
end

% stretch for panel
Z=Z(:,:);

n_endo=size(SIGu,1);

na_over_K=size(V,1);

na=n_endo*na_over_K;

na2=na;

a2=@(x)x;

K=1;

iSIGu=SIGu\eye(n_endo);

dfunc=@()0;

afunc=@(x)x;

if ~isempty(linres)
    
    a2=@(x)linres.a_to_a2tilde(x);
    
    afunc=@(x)linres.a2tilde_to_a(x);
    
    K=linres.K;
    
    dfunc=@()kron(Z*Z.',iSIGu)*linres.d;
    
    na2=size(K,2);
    
end

Sstar=eye(n_endo);

n=n_endo+1;

iV=V\eye(na_over_K);

SIGa=@(SIGu)inv(K.'*(kron(iV+Z*Z.',SIGu\eye(n_endo)))*K); % 5.2.15

abar=SIGa(SIGu)*(K.'*kron(iV,iSIGu)*K*a2(astar)+...
    K.'*(kron(Z,iSIGu)*y-dfunc())); % 5.2.16

% SIGa=@(SIGu)kron(inv(iV+Z*Z.'),SIGu); % 5.2.15
% 
% abar=SIGa(SIGu)*[kron(iV,iSIGu),kron(Z,iSIGu)]*[astar;y(:)]; % 5.2.16

Astar=reshape(astar,n_endo,[]);

Abar=reshape(afunc(abar),n_endo,[]);

Y=reshape(y,n_endo,[]);

% Pseudo-OLS with restrictions
%-----------------------------
ahat=(K.'*(kron(Z*Z.',iSIGu))*K)\(K.'*(kron(Z,iSIGu)*y));

Ahat=reshape(afunc(ahat),n_endo,[]);

% Ahat=Y/Z; % Y*Z.'*inv(Z*Z.')

T=size(Z,2);

tau=T+n;

SIGu_tilde=(Y-Ahat*Z)*(Y-Ahat*Z).'/T;

% I have not managed to understand why the following is usually not
% positive definite
%-------------------------------------------------------------------------
S=T*SIGu_tilde+Sstar+Ahat*(Z*Z.')*Ahat.'+Astar*iV*Astar.'-Abar*(iV+Z*Z.')*Abar.'; % 5.2.17

smpl=@posterior_simulator;

% expand for output!
abar_out=afunc(abar);

    function draws=posterior_simulator(ndraws)
               
        SIGud=SIGu;
        
        for idraw=1:ndraws
            
            CSIGa=chol(SIGa(SIGud),'lower');
            
            % re-expand immediately
            adraw=afunc(abar+CSIGa*randn(na2,1));
            
            di=[adraw;vech(SIGud)];
            
            if idraw==1
                
                draws=di(:,ones(ndraws,1));
                
            else
                
                draws(:,idraw)=di;
                
            end
            
            SIGud=iwishrnd(S,tau);
            
        end
        
    end

end