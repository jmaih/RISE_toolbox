function  [X,retcode]=realize_solution(X,doall)

% this function is about how to treat complex solutions

if nargin<2
    
    doall=false;
    
end

retcode=0;

Xnew=real(X);

if max(abs(X(:)-Xnew(:)))>1e-8
    
    if ~doall
        
        retcode=223; % complex solution
        
    end
    
else
    
    X=Xnew;
    
end

end