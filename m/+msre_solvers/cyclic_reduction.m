function [S,retcode]=cyclic_reduction(Aplus,A0,Aminus,tol,maxiter,debug)
% Aplus*S^2+A0*S+Aminus=0

% Reference: "Efficient cyclic reduction for QBDs with rank structured
% blocks" https://arxiv.org/pdf/1601.00861.pdf

if nargin<6
    
    debug=[];
    
    if nargin<5
        
        maxiter=[];
        
        if nargin<4
            
            tol=[];
            
        end
        
    end
    
end

if isempty(maxiter),maxiter=100;end

if isempty(tol),tol=sqrt(eps);end

if isempty(debug),debug=false;end

Sk=A0;

Ak=Aplus;

Bk=A0;

Ck=Aminus;

done=false;

convfunc=@(X1,X2)max(abs(X1(:)-X2(:)));

iter=0;

n=size(A0,2);

In=eye(n);

retcode=0;

while ~done
    
    iter=iter+1;
    
    iBk=Bk\In;
    
    if any(isnan(iBk(:)))
        
        retcode=22;
        
        done=true;
        
    else
        
        Sk1=Sk-Ak*iBk*Ck; % check
        
        Ak1=-Ak*iBk*Ak; % check
        
        Bk1=Bk-Ak*iBk*Ck-Ck*iBk*Ak; % check
        
        Ck1=-Ck*iBk*Ck; % check
        
        % convergence
        best=convfunc(Ak1,Ak);
        
        if debug
            
            best=max(best,convfunc(Sk,Sk1));
            
            best=max(best,convfunc(Bk,Bk1));
            
            best=max(best,convfunc(Ck,Ck1));
            
        end
        
        % next round
        Sk=Sk1; Ak=Ak1; Bk=Bk1; Ck=Ck1;
        
        done=best<tol;
        
        if debug
            
            fprintf(1,'iter = %0.0f, conv = %0.15f\n',iter,best)
            
        end
        
    end
    
    if iter==maxiter
        
        retcode=21;
        
        done=true;
        
    end
    
end

if retcode
    
    S=[];
    
else
    
    S=-Sk\Aminus;
    
end

end