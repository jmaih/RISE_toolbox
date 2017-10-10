function [abar_out,SIGu,smpl]=normal_prior_with_known_Sigma(Z,SIGu,Va,astar,y,linres)

if nargin < 6
    
    linres=[];
    
end

% stretch for panel
Z=Z(:,:);

n_endo=size(SIGu,1);

na=size(Va,1);

na2=na;

iVa=Va\eye(na);

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

% SIGa=(iVa+kron(Z*Z.',iSIGu))\eye(na); % 5.2.7
% 
% abar=SIGa*(iVa*astar+kron(Z,iSIGu)*y); % 5.2.6

SIGa=(K.'*(iVa+kron(Z*Z.',iSIGu))*K)\eye(na2); % 5.2.7

abar=SIGa*(K.'*iVa*K*a2(astar)+...
    K.'*(kron(Z,iSIGu)*y-dfunc())); % 5.2.6

CSIGa=chol(SIGa,'lower');

% clear iVa iSIGu SIGa

smpl=@posterior_simulator;

% expand for output!
abar_out=afunc(abar);

    function draws=posterior_simulator(ndraws)
        
        for idraw=1:ndraws
            % re-expand immediately
            adraw=afunc(abar+CSIGa*randn(na2,1));
            
            if idraw==1
                
                di=[adraw;vech(SIGu)];
                
                draws=di(:,ones(ndraws,1));
                
            else
                
                draws(1:na,idraw)=adraw;
                
            end
            
        end
        
    end

end