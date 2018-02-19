function pen=penalize_violations2(M,nonlinres,p)

if nargin<3
    
    p=1000;
    
end

ns=size(M,3);

pen=zeros(1,ns);

if ~isempty(nonlinres)
    
    for iii=1:ns
        
        viol=nonlinres(M(:,:,iii));
        
        good_ones=viol<=0;
        
        good=all(good_ones);
        
        if ~good
            
            pen(iii)=uminus(p*sum(viol(~good_ones)));
            
%             retcode(iii)=7;
            
        end
        
    end
    
end
