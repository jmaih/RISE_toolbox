function m=normal(order)

double_factorial=@(n)prod(n:-2:0);

m=zeros(order,1);

for o=1:order
    
    if rem(o,2)==0
        
        m(o)=double_factorial(o-1);
        
    end
    
end

end
