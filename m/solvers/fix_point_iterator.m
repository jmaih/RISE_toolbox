%  fix_point_iterator solves the fix point of a function
% 
%  ::
% 
%    [T,itercode,retcode]=fix_point_iterator(iterate_func,T0,options,varargin)
% 
%  Args:
% 
%     - iterate_func : [func_handle]: a function handle which returns 2
%       elements: The improved value of T and F0, the value of the function
%       whose discrepancy is to be minmized
% 
%     - T0 : is the initial guess for the solution
% 
%     - options : [struct] with fields
%       - fix_point_explosion_limit : [positive scalar |{1e+6}] : maximum
%         divergence of F0
%       - fix_point_TolFun : [positive scalar |{sqrt(eps)}] : tolerance
%         criterion
%       - fix_point_maxiter : [positive scalar |{1000}] : maximum number of
%         iterations
%       - fix_point_verbose : [true|{false}] : show iterations or not
%       - fix_point_valid_func : [function_handle|{@utils.error.valid}] :
%         function that checks the validity (nans, complex, inf, etc.) of
%         the calculations. The default function checks for complex numbers.
%         If you want to allow for complex numbers, use instead
%         @utils.error.validComplex
% 
%     - varargin
% 
%  Returns:
%     :
%     - T : final solution
%     - itercode : final number of iterations
%     - retcode : return code
%       - 0 : successful
%       - 21 : maximum number of iterations reached
%       - 22 : nan or inf in F0 or T
%       - 23 : divergence
% 
%  Note:
% 
%  Example:
% 
%     See also:
%