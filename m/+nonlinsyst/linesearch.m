function [dx,fvalnew,fnormnew,flag,iter]=linesearch(f,x,dx,fnorm,maxsteps,varargin)

fnormold = inf;

iter=0;

while iter<maxsteps
    
    iter=iter+1;
    
    fvalnew = f(x+dx,varargin{:});
    
    fnormnew = nonlinsyst.norm(fvalnew);
    
    if fnormnew<fnorm
        
        flag=0;
        
        return
        
    end
    
    if fnormold<fnormnew
        
        fvalnew=fvalold;
        
        fnormnew=fnormold;
        
        dx=dx*2;
        
        flag=2;
        
        return
        
    end
    
    fvalold  = fvalnew;
    
    fnormold = fnormnew;
    
    dx = dx/2;
    
end

flag=1;

end