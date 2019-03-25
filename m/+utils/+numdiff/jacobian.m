%--- help for sym/jacobian ---
%
% JACOBIAN Jacobian matrix.
%    JACOBIAN(f,x) computes the Jacobian of the scalar or vector f
%    with respect to the vector x. The (i,j)-th entry of the result
%    is df(i)/dx(j). Note that when f is scalar, the Jacobian of f
%    is the gradient of f. Also, note that scalar x is allowed,
%    although this is just DIFF(f,x).
% 
%    Example:
%        syms x y z u v; jacobian([x*y*z; y; x+z],[x y z])
%        returns  [y*z, x*z, x*y; 0, 1, 0; 1, 0, 1]
% 
%        jacobian(u*exp(v),[u;v])
%        returns  [exp(v), u*exp(v)]
% 
%    See also SYM/CURL, SYM/DIFF, SYM/DIVERGENCE, SYM/GRADIENT, SYM/HESSIAN,
%    SYM/POTENTIAL, CURL, DIVERGENCE, HESSIAN, LAPLACIAN, VECTORPOTENTIAL,
%    SUBS.
%