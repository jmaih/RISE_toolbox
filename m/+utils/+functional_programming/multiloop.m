function [b,out]=multiloop(b,start,pos,start_init,dofunc,out,varargin)
% this function is designed to replace multi for loops
% b : vector gathering the indices of the for the different for loops
% start: index where the current for loop will start
% pos: position of the current for loop
% start_init: true or false, start the next for loop from 1 or from the
% same index as the previous one
% dofunc : at the end of one pass, i.e. b is completed, it is evaluated.
% This is a function handle
% out : input and output of the dofunc function

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