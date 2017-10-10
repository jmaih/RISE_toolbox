function y=simulate(y0,xdet,A,B,shocks)

[nshocks,nperiods]=size(shocks);

nvars=nshocks;

y=zeros(nvars,nperiods);

y0t=y0;

ndet=size(xdet,2);

for t=1:nperiods
    
    tdet=min(t,ndet);
    
    y(:,t)=A*flip_initial(y0t)+B*shocks(:,t);
    
    y0t=[y0t(:,2:end),y(:,t)];
    
end

    function xx=flip_initial(x)
        
        xx=fliplr(x); % xx=x(end:-1:1)
        
        xx=xx(:);
        
        if ~isempty(xdet)
            
            xx=[xdet(:,tdet);xx];
            
        end
        
    end

end
