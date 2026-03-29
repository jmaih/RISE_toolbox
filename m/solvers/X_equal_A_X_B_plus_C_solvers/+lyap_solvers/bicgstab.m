% BICGSTAB   BiConjugate Gradients Stabilized Method.
%    X = BICGSTAB(A,B) attempts to solve the system of linear equations
%    A*X=B for X. The N-by-N coefficient matrix A must be square and the
%    right hand side column vector B must have length N.
% 
%    X = BICGSTAB(AFUN,B) accepts a function handle AFUN instead of the
%    matrix A. AFUN(X) accepts a vector input X and returns the
%    matrix-vector product A*X. In all of the following syntaxes, you can
%    replace A by AFUN.
% 
%    X = BICGSTAB(A,B,TOL) specifies the tolerance of the method. If TOL is
%    [] then BICGSTAB uses the default, 1e-6.
% 
%    X = BICGSTAB(A,B,TOL,MAXIT) specifies the maximum number of iterations.
%    If MAXIT is [] then BICGSTAB uses the default, min(N,20).
% 
%    X = BICGSTAB(A,B,TOL,MAXIT,M) and X = BICGSTAB(A,B,TOL,MAXIT,M1,M2) use
%    preconditioner M or M=M1*M2 and effectively solve the system
%    A*inv(M)*Y = B for Y, where Y = M*X. If M is [] then a preconditioner
%    is not applied. M may be a function handle returning M\X.
% 
%    X = BICGSTAB(A,B,TOL,MAXIT,M1,M2,X0) specifies the initial guess.  If
%    X0 is [] then BICGSTAB uses the default, an all zero vector.
% 
%    [X,FLAG] = BICGSTAB(A,B,...) also returns a convergence FLAG:
%     0 BICGSTAB converged to the desired tolerance TOL within MAXIT iterations.
%     1 BICGSTAB iterated MAXIT times but did not converge.
%     2 preconditioner M was ill-conditioned.
%     3 BICGSTAB stagnated (two consecutive iterates were the same).
%     4 one of the scalar quantities calculated during BICGSTAB became
%       too small or too large to continue computing.
% 
%    [X,FLAG,RELRES] = BICGSTAB(A,B,...) also returns the relative residual
%    NORM(B-A*X)/NORM(B). If FLAG is 0, then RELRES <= TOL.
% 
%    [X,FLAG,RELRES,ITER] = BICGSTAB(A,B,...) also returns the iteration
%    number at which X was computed: 0 <= ITER <= MAXIT. ITER may be an
%    integer + 0.5, indicating convergence half way through an iteration.
% 
%    [X,FLAG,RELRES,ITER,RESVEC] = BICGSTAB(A,B,...) also returns a vector
%    of the residual norms at each half iteration, including NORM(B-A*X0).
% 
%    Example:
%       n = 21; A = gallery('wilk',n);  b = sum(A,2);
%       tol = 1e-12;  maxit = 15; M = diag([10:-1:1 1 1:10]);
%       x = bicgstab(A,b,tol,maxit,M);
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
%    as inputs to BICGSTAB:
%       x1 = bicgstab(@(x)afun(x,n),b,tol,maxit,@(x)mfun(x,n));
% 
%    Class support for inputs A,B,M1,M2,X0 and the output of AFUN:
%       float: double
% 
%    See also BICG, BICGSTABL, CGS, GMRES, LSQR, MINRES, PCG, QMR, SYMMLQ,
%    TFQMR, ILU, FUNCTION_HANDLE.
%
%    Documentation for bicgstab
%       doc bicgstab
%
%    Other uses of bicgstab
%
%       codistributed/bicgstab    gpuArray/bicgstab
%