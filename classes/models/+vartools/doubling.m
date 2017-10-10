function [P]=doubling(A,B,C,fix_point_verbose)
% doubling_solve solves the linear equation X=A*X*B+C
%
% Syntax
% -------
% ::
%   [P,retcode]=doubling(A,B,C)
%

if nargin < 4
    
    fix_point_verbose=[];
    
    if nargin<3
        
        error([mfilename,':: at least 3 arguments should be provided'])
        
    end
    
end

if isempty(fix_point_verbose)
    
    fix_point_verbose=false;
    
end

fix_point_TolFun=sqrt(eps);

fix_point_maxiter=1000;

fix_point_explosion_limit=1e+12;

valid=@(x)isreal(x) && ~any(isnan(x(:))) && ~any(isinf(x(:)));

if isempty(B),B=A';end

if isempty(A),A=B';end


P0=C;

symmetric=isequal(A,B');

Gl=A;

if ~symmetric
    
    Gr=B;
    
end

iter=0;

conv_F=0.5*fix_point_explosion_limit;

conv_T=conv_F;

F0=rand;

while max(conv_T,conv_F)>fix_point_TolFun && ...
        iter<fix_point_maxiter && ...
        conv_F<fix_point_explosion_limit && ...
        valid(F0)
    
    iter=iter+1;
    
    [P,F0] = iterator(P0);
    
    conv_T=max(abs(P(:)-P0(:)));
    
    conv_F=max(abs(F0(:)));
    
    if fix_point_verbose
        
        fprintf(1,'iter # %0.0f : conv(x)=%0.4f, conv(F(x))=%0.4f\n',iter,full(conv_T),conv_F);
        
    end
    
    P0=P;
    
    if ~valid(P)
        
        break
        
    end
    
end

    function [P,F0]=iterator(P0)
        
        if symmetric
            
            P=P0+Gl*P0*Gl';
            
        else
            
            P=P0+Gl*P0*Gr;
            
            Gr=Gr*Gr;
            
        end
        
        Gl=Gl*Gl;
        
        F0=P-P0;
        
    end

end