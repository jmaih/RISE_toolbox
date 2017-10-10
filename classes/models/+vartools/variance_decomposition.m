function [Vinfi,Vi,info]=variance_decomposition(A,B,nz,max_periods)

if nargin<4
    
    max_periods=[];
    
end

if isempty(max_periods),max_periods=400; end

% collect the state matrices
endo_nbr=size(A,1);

k=max(100,max(max_periods));

exo_nbr=endo_nbr;

[T,R]=vartools.companion(A,B,nz);

[Vinfi,Vi]=theoretical_vardec_engine(T,R);

Vinfi=Vinfi(1:endo_nbr,:,:);

Vi=Vi(1:endo_nbr,:,:);

info={'endogenous','shocks','horizon'};

    function [Vinfi,Vi]=theoretical_vardec_engine(T,R)
        
        grand_endo_nbr=size(R,1);
        
        R=reshape(R,grand_endo_nbr,exo_nbr);
        
        V=zeros(grand_endo_nbr);
        
        Vi=zeros(grand_endo_nbr,exo_nbr,k);
        
        Vii=zeros(grand_endo_nbr,grand_endo_nbr,exo_nbr);
        
        RR=R*R';
        
        for ii=1:k
            
            V=T*V*T'+RR;
            
            [Vi(:,:,ii),Vii]=decompose_variance(Vi(:,:,ii),V,Vii);
            
        end
        
        Vinfi=zeros(grand_endo_nbr,exo_nbr);
        
        [Vinf]=vartools.solve_lyapunov_equation(T,RR);
        
        if any(~isfinite(Vinf(:)))
            
            error('Variance could not be solved')
            
        end
        
        % deal with zero variances
        Vinfi=decompose_variance(Vinfi,Vinf);
        
        function [V,Vkk]=decompose_variance(V,total_variance,Vkk)
            
            total_variance=diag(total_variance);
            
            total_variance(total_variance<1e-12)=1;
            
            for iexo=1:exo_nbr
                
                Ri=zeros(size(R));
                
                locs=iexo:exo_nbr:exo_nbr;
                
                Ri(:,locs)=R(:,locs);
                
                RRi=Ri*Ri';
                
                if nargin<3
                    
                    [V00]=vartools.solve_lyapunov_equation(T,RRi);
                    
                    if any(~isfinite(V00(:)))
                        
                        error('Variance could not be solved')
                        
                    end
                    
                    Vkk=V00;
                    
                else
                    
                    Vkk(:,:,iexo)=T*Vkk(:,:,iexo)*T'+RRi;
                    
                    V00=Vkk(:,:,iexo);
                    
                end
                
                V(:,iexo)=diag(V00)./total_variance;
                
            end
            
        end
        
    end

end
