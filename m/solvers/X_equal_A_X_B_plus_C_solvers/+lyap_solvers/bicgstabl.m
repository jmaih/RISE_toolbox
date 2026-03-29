% BICGSTABL   BiConjugate Gradients Stabilized(l) Method.
%    X = BICGSTABL(A,B) attempts to solve the system of linear equations
%    A*X=B for X. The N-by-N coefficient matrix A must be square and the
%    right hand side column vector B must have length N.
% 
%    X = BICGSTABL(AFUN,B) accepts a function handle AFUN instead of the
%    matrix A. AFUN(X) accepts a vector input X and returns the
%    matrix-vector product A*X. In all of the following syntaxes, you can
%    replace A by AFUN.
% 
%    X = BICGSTABL(A,B,TOL) specifies the tolerance of the method. If TOL is
%    [] then BICGSTABL uses the default, 1e-6.
% 
%    X = BICGSTABL(A,B,TOL,MAXIT) specifies the maximum number of iterations.
%    If MAXIT is [] then BICGSTABL uses the default, min(N,20).
% 
%    X = BICGSTABL(A,B,TOL,MAXIT,M) and X = BICGSTABL(A,B,TOL,MAXIT,M1,M2)
%    use preconditioner M or M=M1*M2 and effectively solve the system
%    A*inv(M)*Y = B for Y, where Y = M*X. If M is [] then a preconditioner
%    is not applied. M may be a function handle returning M\X.
% 
%    X = BICGSTABL(A,B,TOL,MAXIT,M1,M2,X0) specifies the initial guess.  If
%    X0 is [] then BICGSTABL uses the default, an all zero vector.
% 
%    [X,FLAG] = BICGSTABL(A,B,...) also returns a convergence FLAG:
%     0 BICGSTABL converged to the desired tolerance TOL within MAXIT iterations.
%     1 BICGSTABL iterated MAXIT times but did not converge.
%     2 preconditioner M was ill-conditioned.
%     3 BICGSTABL stagnated (two consecutive iterates were the same).
%     4 one of the scalar quantities calculated during BICGSTABL became
%       too small or too large to continue computing.
% 
%    [X,FLAG,RELRES] = BICGSTABL(A,B,...) also returns the relative residual
%    NORM(B-A*X)/NORM(B). If FLAG is 0, then RELRES <= TOL.
% 
%    [X,FLAG,RELRES,ITER] = BICGSTABL(A,B,...) also returns the iteration
%    number at which X was computed: 0 <= ITER <= MAXIT. ITER may be k/4 where
%    k is some integer, indicating convergence at a given quarter iteration.
% 
%    [X,FLAG,RELRES,ITER,RESVEC] = BICGSTABL(A,B,...) also returns a vector
%    of the residual norms at each quarter iteration, including NORM(B-A*X0).
% 
%    Example:
%       n = 21; A = gallery('wilk',n);  b = sum(A,2);
%       tol = 1e-12;  maxit = 15; M = diag([10:-1:1 1 1:10]);
%       x = bicgstabl(A,b,tol,maxit,M);
%    Or, use this matrix-vector product function
%       %-----------------------------------------------------------------%
%       function y = afun(x,n)
%       y = [0; x(1:n-1)] + [((n-1)/2:-1:0)'; (1:(n-1)/2)'].*x+[x(2:n); 0];
%       %-----------------------------------------------------------------%
%    and this preconditioner backsolve function
%       %------------------------------------------%
%       function y = mfun(r,n)
%       y = r ./ [((n-1)/2:-1:1)'; 1; (1:(n-1)/2)'];
%       %------------------------------------------%
%    as inputs to BICGSTABL:
%       x1 = bicgstabl(@(x)afun(x,n),b,tol,maxit,@(x)mfun(x,n));
% 
%    Class support for inputs A,B,M1,M2,X0 and the output of AFUN:
%       float: double
% 
%    See also BICG, BICGSTAB, CGS, GMRES, LSQR, MINRES, PCG, QMR, SYMMLQ,
%    TFQMR, ILU, FUNCTION_HANDLE.
%
%    Documentation for bicgstabl
%       doc bicgstabl
%
%    Other uses of bicgstabl
%
%       codistributed/bicgstabl    gpuArray/bicgstabl
%