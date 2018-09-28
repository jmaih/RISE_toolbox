function M=set_up_normal_moments(order,nshocks)

% https://www.dsprelated.com/freebooks/sasp/Gaussian_Moments.html

double_factorial=@(n)prod(n:-2:0);

M=cell(1,order);

for oo=1:order
    
    M{oo}=set_up_one(oo);
    
end

    function m=set_up_one(o)
        
        if o==1
            
            m=sparse(nshocks,1);
            
        else
            
            siz=nshocks*ones(1,o);
            
            m=zeros(siz);
            
            if rem(o,2)==0
                
                e=repmat({(1:nshocks).'},1,o);
                
                IND=sub2ind(siz,e{:});
                
                m(IND)=double_factorial(o-1);
                
            end
            
        end
        
    end

end
