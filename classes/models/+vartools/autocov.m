function [Vi,info]=autocov(A,B,nz,max_periods)

if nargin<4
    
    max_periods=[];
    
end

if isempty(max_periods),max_periods=5; end

% collect the state matrices
endo_nbr=size(A,1);

k=max_periods;

exo_nbr=endo_nbr;

[T,R]=vartools.companion(A,B,nz);

Vi=autocov_engine(T,R);

Vi=Vi(1:endo_nbr,1:endo_nbr,:);

info={'endogenous','endogenous','horizon'};

    function Vi=autocov_engine(T,R)
        
        grand_endo_nbr=size(R,1);
        
        R=reshape(R,grand_endo_nbr,exo_nbr);
        
        Vi=zeros(grand_endo_nbr,grand_endo_nbr,k+1);
                
        RR=R*R';
        
        [V]=vartools.solve_lyapunov_equation(T,RR);
        
        if any(~isfinite(V(:)))
            
            error('Variance could not be solved')
            
        end
        
        for ii=1:k+1
            
            Vi(:,:,ii)=V;
            
            V=T*V;
            
        end
                
    end

end