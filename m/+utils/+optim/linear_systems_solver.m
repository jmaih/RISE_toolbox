function [X,retcode]= linear_systems_solver(A,b,x0,options)
% linear_systems_solver -- solver for linear sytems
%
% ::
%
%
%   [X,retcode]= linear_systems_solver(A,b,x0,options)
%
% Args:
%
%    - **A** [matrix|function_handle]: such that lhs=A*x or lhs=A(x)
%
%    - **b** [vector]: right hand side of the system
%
%    - **x0** [vector]: initial guess of the solution
%
%    - **options** [struct]: rise options. The most relevant are
%      - **fix_point_TolFun** [numeric]: tolerance level
%      - **fix_point_maxiter** [integer]: maximum number of iterations
%      - **fix_point_verbose** [true|false]: display iterations (default
%      solver only)
%      - **solve_linsyst_user_algo** [{[]}|char|function_handle|cell]:
%      user-defined linear solver. The solver must take the same inputs and
%      outputs as Matlab solvers: tfqmr, bicg, bicgstab, bicgstabl, cgs,
%      minres, pcg, qmr. When empty, the default solver
%      transpose_free_quasi_minimum_residual is applied
%
% Returns:
%    :
%
%    - **X** [vector]: solution
%
%    - **retcode** [numeric]: if retcode=0, then a solution is found. Else
%    there is a problem whose meaning can be found by running
%    decipher(retcode)
%
% Note:
%
%    - Matlab functions tfqmr, bicg, bicgstab, bicgstabl, cgs, minres, pcg,
%    qmr can be readily used. In order to use gmres, the user has to define a
%    function handle since gmres does not have the same inputs as the others.
%    E.g. one can define
%    gmres_=@(A,b,tol,maxit,M1,M2,x0,varargin)gmres(A,b,restart,tol,maxit,M1,M2,x0,varargin{:})
%    where restart has been predefined.
%
%    - Inputs M1 and M2 are always empty i.e. there is no preconditioning
%
% Example:
%
%    See also:

if ~isfield(options,'solve_linsyst_user_algo')
    
    options.solve_linsyst_user_algo=[];
    
end

if isempty(options.solve_linsyst_user_algo)
    [X,retcode]=transpose_free_quasi_minimum_residual(A,b,... % right hand side
        x0,... % initial guess
        options.fix_point_TolFun,... % tolerance level
        options.fix_point_maxiter,... % maximum number of iterations
        options.fix_point_verbose);
else
    % [x,flag,relres,iter,resvec] = tfqmr(A,b,tol,maxit,M1,M2,x0,varargin)
    % [x,flag,relres,iter,resvec] = bicg(A,b,tol,maxit,M1,M2,x0,varargin)
    % [x,flag,relres,iter,resvec] = bicgstab(A,b,tol,maxit,M1,M2,x0,varargin)
    % [x,flag,relres,iter,resvec] = bicgstabl(A,b,tol,maxit,M1,M2,x0,varargin)
    % [x,flag,relres,iter,resvec] = cgs(A,b,tol,maxit,M1,M2,x0,varargin)
    % [x,flag,relres,iter,resvec] = gmres(A,b,restart,tol,maxit,M1,M2,x,varargin)
    % [x,flag,relres,iter,resvec,resveccg] = minres(A,b,tol,maxit,M1,M2,x0,varargin)
    % [x,flag,relres,iter,resvec] = pcg(A,b,tol,maxit,M1,M2,x0,varargin)
    % [x,flag,relres,iter,resvec] = qmr(A,b,tol,maxit,M1,M2,x0,varargin)
    M1=[];
    M2=[];
    [linsolver,vargs]=utils.code.user_function_to_rise_function(...
        options.solve_linsyst_user_algo);
    [X,flag,relres,iter,resvec] = linsolver(A,b,...
        options.fix_point_TolFun,...
        options.fix_point_maxiter,M1,M2,x0,vargs{:}); %#ok<ASGLU>
    if flag==0
        retcode=0;
    elseif flag==1
        % maximum number of iterations reached
        retcode=201;
    else
        % Nans in solution or no solution
        retcode=202;
    end
end

end