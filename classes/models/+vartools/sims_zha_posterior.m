function [abar_out,SIGu,smpl]=sims_zha_posterior(Z,SIGu,omega,y,linres)

if nargin<5
    
    linres=[];
    
end

% stretch for panel
Z=Z(:,:);

n_endo=size(SIGu,1);

Y=reshape(y,n_endo,[]);

na=n_endo*size(Z,1);% na=size(Va,1);

na2=na;

% iVa=Va\eye(na);

% a2=@(x)x;

K=1;

iSIGu=eye(n_endo);

dfunc=@()0;

afunc=@(x)x;

if ~isempty(linres)
    
%     a2=@(x)linres.a_to_a2tilde(x);
    
    afunc=@(x)linres.a2tilde_to_a(x);
    
    K=linres.K;
    
    dfunc=@()kron(Z*Z.',iSIGu)*linres.d;
    
    na2=size(K,2);
    
end

[~,abar]=calculate_SIGa_a(eye(n_endo));

% re-expand immediately
abar_out=afunc(abar);

B=reshape(abar_out,n_endo,[]);

resids=Y-B*Z;

T=size(Y,2);

Tminus=0; % training sample

tau=n_endo*(omega+1)+Tminus+T;

S=resids*resids.';

SIGu=S/T;

[SIGa]=calculate_SIGa_a(SIGu);

smpl=@posterior_simulator;
            

    function [SIGa,abar]=calculate_SIGa_a(SIGu)
        
        iSIGu=SIGu\eye(n_endo);
                
        SIGa=(K.'*(kron(Z*Z.',iSIGu))*K)\eye(na2); % 5.2.24
        
        if nargout>1
            
            abar=SIGa*(K.'*(kron(Z,iSIGu)*y-dfunc())); % 5.2.23
            
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
            
            SIGud=iwishrnd(S,tau);
            
            SIGad=calculate_SIGa_a(SIGud);
            
        end
        
    end

end