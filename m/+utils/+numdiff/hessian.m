%--- help for sym/hessian ---
%
% HESSIAN Hessian matrix.
%    HESSIAN(f,x) computes the Hessian of the scalar f with respect
%    to the vector x. The (i,j)-th entry of the resulting matrix is
%    (d^2f/(dx(i)dx(j)). Note that scalar x is allowed although this
%    is just DIFF(f,x,2).
% 
%    Example:
%        syms x y z; hessian(x*y + 2*z*x, [x, y, z])
%        returns  [0, 1, 2; 1, 0, 0; 2, 0, 0]
% 
%    See also SYM/CURL, SYM/DIFF, SYM/GRADIENT, SYM/JACOBIAN, SYM/DIVERGENCE,
%    SYM/POTENTIAL, CURL, DIVERGENCE, HESSIAN, JACOBIAN, LAPLACIAN,
%    VECTORPOTENTIAL, SUBS.
%
%    Other functions named hessian
%
%       generic/hessian
%