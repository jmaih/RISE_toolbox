function sol=polynsolve(polyset,ord,varnames,tol)
% POLYNSOLVE - attempt exact solution of multivariate polynomial system
%
% usage: sol=polynsolve(polyset)
%
% INPUTS: polyset, ord, varnames, tol
%  polyset is a cell array of polynomials in string or coefficient form
%    that is acceptable input for petschel.groebner.m
%  ord (optional) is the preferred ordering
%  varnames (optional) is the list of variable names if not {'x1','x2',...}
%  tol (optional) is the default zero tolerance
%    (see petschel.groebner.m for further details on ord, varnames, tol)
%
% OUTPUTS: sol
%  sol is an array containing the solutions to {polyset{:}=0}:
%    sol(i,j) is the value of xj in the i'th solution
%    If there are infinitely many solutions, sol(i,j)=NaN
%
% ALGORITHM:
%  Uses petschel.groebner bases (lex order).  If any one-variable polynomials
%  result, solve them and substitute back into the equations.
%  If "1" results, the system is not solvable.  If any multivariate polys
%  remain, those variables have an infinite number of solutions.
%
% KNOWN BUGS
%  See petschel.groebner.m for details
%
% SEE ALSO:
%  petschel.groebner, petschel.poly2str, petschel.str2poly

% Author: Ben Petschel 23/6/2009
%
% Change history:
%  23/6/2009 - first release
%  30/3/2010 - change poly representation from n-dim to rectangular array

if nargin<4
    
    tol = 0;
    
    if nargin<3
        
        varnames = {};
        
    end
    
end

if (nargin<2) || isempty(ord)
    
    ord = 'lex';
    
end

if (numel(polyset)>0) && ischar(polyset{1})
    
    polyset = petschel.str2poly(polyset,varnames);
    
end

gbasis=petschel.groebner(polyset,ord,varnames,tol);

% if petschel.groebner is all linear terms, it is a solution, otherwise find
% 1-variable polynomials and solve them, substituting the solutions into
% the other equations
sol = [];

i=1;

keepgoing = true;

while keepgoing && (i<=numel(gbasis))
    
    % search for 1-variable polynomials
    d = size(gbasis{i},2)-1;
    
    if numel(sol)<d
        % sol will always be a row vector, until possibly the final step
        % make sure sol has enough possible variables
        sol = [sol,nan(1,d-numel(sol))];
        
    end
    
    [tf,n,Q]=ispoly1(gbasis{i});
    
    if tf
        % attempt solution
        if length(Q)==1
            % have equation 1==0, so no solution
            sol = [];
            
            keepgoing = false;
            
        elseif length(Q)==2
            % Q=[a;b] so have equation a*xn+b=0 (a~=0) and solution is x=-b/a
            
            sol(n) = -Q(2)/Q(1);
            
        else
            % have polynomial in xn, so solve and substitute solution
            r = roots(Q);
            
            gbasisrecur = gbasis;
            
            sol = [];
            
            for j=1:length(r)
                % replace equation i with linear term (xn-r(j)==0)
                
                gbasisrecur{i} = linearterm(r(j),n);
                
                sol = [sol; petschel.polynsolve(gbasisrecur,ord,varnames,tol)];
                
                keepgoing = false;
                
            end
            
        end
        
    end
    
    i = i+1;
    
end


end % main function polynsolve(...)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tf,n,Q]=ispoly1(P)
% returns tf=true if P is a 1-variable polynomial
% n is the number of the variable the polynomial is in
% Q is a column vector of the polynomial coefficients

if isempty(P) || all(P(:)==0)
    
    tf=true;
    
    n=1;
    
    Q=0;
    
elseif size(P,2)==1
    
    tf=true;
    
    n=1;
    
    Q=P;
    
else
    
    ind = find(any(P(:,2:end)>0,1));
    
    if numel(ind)==1
        % is a polynomial in 1 variable if exponents of only 1 variable are >0
        tf=true;
        
        n=ind;
        
        Q=[]; % collect coefficients of the polynomial
        
        Q(P(:,ind+1)+1)=P(:,1);
        
        Q=Q(end:-1:1); % reverse order to be consistent with ROOTS
        
    else
        
        tf=false;
        
        n=0;
        
        Q=[];
        
    end
    
end

end % helper function ispoly1(...)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Q=linearterm(a,n)
% returns coefficient array of a linear term (xn-a)

Q=zeros(2,n+1);

Q(1,1)=-a;

Q(2,1)=1;

Q(2,end)=1;

end % helper function linearterm(...)
