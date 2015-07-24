function [X,retcode]= linear_systems_solver(A,b,x0,options)
% H1 line
%
% Syntax
% -------
% ::
%
% Inputs
% -------
%
% Outputs
% --------
%
% More About
% ------------
%
% Examples
% ---------
%
% See also: 

if isempty(options.solve_linsyst_user_algo)
    [X,retcode]=transpose_free_quasi_minimum_residual(A,b,... % right hand side
        x0,... % initial guess
        options.fix_point_TolFun,... % tolerance level
        options.fix_point_maxiter,... % maximum number of iterations
        options.fix_point_verbose);
else
    [linsolver,vargs]=utils.code.user_function_to_rise_function(...
        options.solve_linsyst_user_algo);
    [X,flag,relres,iter,resvec] = linsolver(A,b,...
        options.fix_point_TolFun,...
        options.fix_point_maxiter,vargs{:}); %#ok<ASGLU>
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