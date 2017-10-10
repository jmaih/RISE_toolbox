function [abar_out,SIGu,smpl]=independent_normal_wishart_posterior(Z,SIGu,Va,...
    astar,y,linres)

if nargin<6
    
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

Sstar=eye(n_endo);

n=n_endo+1;

Y=reshape(y,n_endo,[]);

T=size(Z,2);

[SIGa,abar]=calculate_SIGa_a(SIGu);

tau=T+n;

smpl=@posterior_simulator;

% re-expand immediately
abar_out=afunc(abar);
            

    function [SIGa,abar]=calculate_SIGa_a(SIGu)
        
        iSIGu=SIGu\eye(n_endo);
        
        SIGa=(K.'*(iVa+kron(Z*Z.',iSIGu))*K)\eye(na2); % 5.2.24
        
        if nargout>1
            
            abar=SIGa*(K.'*iVa*K*a2(astar)+...
                K.'*(kron(Z,iSIGu)*y-dfunc())); % 5.2.23
        end

    end

    function draws=posterior_simulator(ndraws)
               
        SIGud=SIGu;
        
        SIGad=SIGa;
        
        for idraw=1:ndraws
            
            CSIGa=chol(SIGad,'lower');
            
            % re-expand immediately
            adraw=afunc(abar+CSIGa*randn(na2,1));
            
            di=[adraw;vech(SIGud)];
            
            if idraw==1
                
                draws=di(:,ones(ndraws,1));
                
            else
                
                draws(:,idraw)=di;
                
            end
            
            S=calculate_S();
            
            SIGud=iwishrnd(S,tau);
            
            SIGad=calculate_SIGa_a(SIGud);
            
        end
        
        function S=calculate_S()
            
            Ad=reshape(adraw,n_endo,[]);
            
            yaz=Y-Ad*Z;
            
            S=Sstar+(yaz*yaz.');
            
        end
        
    end

end