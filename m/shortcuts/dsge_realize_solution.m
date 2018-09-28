function  [X,retcode]=dsge_realize_solution(X,doall)

% this function is about how to treat complex solutions

if nargin<2
    
    doall=false;
    
end

retcode=0;

if doall
    
    return
    
end

% Xr=real(X); Xi=imag(X); Xi(abs(Xi)<1e-8)=0; X=Xr+Xi*1i;

Xold=X;

X=real(X);

if max(abs(X(:)-Xold(:)))>1e-8
    
    retcode=223; % complex solution
    
end

end