function [T,retcode,T2]=refine(T1,A0,Aplus,Q,opts)
% INTERNAL FUNCTION
%

% T1 [m x m2 x h array]: previous solution when sig=0
% A0 [1 x h cell] : coefficient on contemporaneous terms
% Aplus [1 x h cell] : coefficient on future terms when  sig=0
% Q [h x h matrix] : Transition matrix at the steady state

h=size(Q,1);

m=size(T1,1);

PHI=zeros(m,m,h);

Im=eye(m);

Abarplus=cell(h);

for ii=1:h
    
    pbar_i=1-Q(ii,ii);
    
    A0AiTi=A0{ii}+(1-pbar_i)*Aplus{ii}*T1(:,:,ii);
    
    iA0AiTi=A0AiTi\Im;
    
    PHI(:,:,ii)=pbar_i*Aplus{ii}*T1(:,:,ii)^2;
    
    for jj=1:h
        
        if ii==jj
            
            continue
            
        end
        
        AplusTilde_ij=Aplus{ii}*Q(ii,jj);
        
        PHI(:,:,ii)=PHI(:,:,ii)-AplusTilde_ij*T1(:,:,jj)*T1(:,:,ii);
        
        Abarplus{ii,jj}=iA0AiTi*AplusTilde_ij;
        
    end
    
    PHI(:,:,ii)=iA0AiTi*PHI(:,:,ii);
    
    Abarplus{ii,ii}=iA0AiTi*(1-pbar_i)*Aplus{ii};
    
end

T20=zeros(m,m,h);

[T2,~,retcode]=fix_point_iterator(@func_iterate,T20(:),opts);

T=T1+T2;

    function [T2,FVAL]=func_iterate(T20)
        
        T20=reshape(T20,m,m,h);
        
        T2=PHI;
        
        for iii=1:h
            
            T2(:,:,iii)=T2(:,:,iii)-Abarplus{iii,iii}*...
                (T20(:,:,iii)*T1(:,:,iii)+T20(:,:,iii)^2);
            
            for jjj=1:h
                
                if jjj==iii
                    
                    continue
                    
                end
                
                T2(:,:,iii)=T2(:,:,iii)-Abarplus{iii,jjj}*...
                    (T1(:,:,jjj)*T20(:,:,iii)+T2(:,:,jjj)*T1(:,:,iii)+...
                    T20(:,:,jjj)*T20(:,:,iii));                    
                
            end
            
        end
        
        FVAL=max(max(abs(T2-T20)));
        
    end

end

% function [T,T2]=refine_maih_waggoner_solution(T1,A0,Aplus,Q)
% % INTERNAL FUNCTION
% %
% 
% % T1 [m x m2 x h array]: previous solution when sig=0
% % A0 [1 x h cell] : coefficient on contemporaneous terms
% % Aplus [1 x h cell] : coefficient on future terms when  sig=0
% % Q [h x h matrix] : Transition matrix at the steady state
% 
% m=size(A0{1},1);
% 
% % Im=eye(m);
% 
% h=numel(A0);
% 
% [pijAiTi,pAplus01]=create_definitions_and_precondition();
% 
% afun=@iterate_func1;
% 
% T2 = tfqmr(afun,pijAiTi(:));
% 
% T2=reshape(T2,[m,m,h]);
% 
% T=T1+T2;
% 
%     function pijAijT2=iterate_func1(T2)
%         
%         T2=reshape(T2,[m,m,h]);
%         
%         pijAijT2=zeros(m,m,h);
%         
%         for r0=1:h
%             
%             for r1=1:h
%                 
%                 if r1==r0
%                     
%                     continue
%                     
%                 end
%                 
%                 pijAijT2(:,:,r0)=pijAijT2(:,:,r0)+pAplus01{r0,r1}*T2(:,:,r1);
%                 
%             end
%             
%         end
%         
%         pijAijT2=pijAijT2(:);
%         
%     end
% 
%     function [pijAiTi,pAplus01]=create_definitions_and_precondition()
%         
%         pijAiTi=zeros(m,m,h);
%         
%         pAplus01=cell(h);
%         
%         for r0=1:h
%                         
%             pijAiTi(:,:,r0)=(1-Q(r0,r0))*Aplus{r0}*T1(:,:,r0);
%             
%             for r1=1:h
%                 
%                 pAplus01{r0,r1}=Q(r0,r1)*Aplus{r0};
%                 
%             end
%             
%         end
%         
%     end
% 
% end