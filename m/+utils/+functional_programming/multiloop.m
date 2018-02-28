function [b,out]=multiloop(b,start,pos,start_init,dofunc,out,varargin)

n=size(out,2);

for ii=start:n
    
    b(pos)=ii;
    
    if pos==length(b)
    
        out=dofunc(out,b,varargin{:});
    
    else
        
        newstart=ii;
        
        if start_init
            
            newstart=1;
            
        end
        
        [b,out]=utils.functional_programming.multiloop(b,newstart,pos+1,...
            start_init,dofunc,out,varargin{:});
        
    end
    
end

end