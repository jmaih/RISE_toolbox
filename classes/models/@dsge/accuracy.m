%--- help for dsge/accuracy ---
%
%  accuracy : computes the Euler errors of dsge models
% 
%  ::
% 
%    Euler_errors=accuracy(m)
%    Euler_errors=accuracy(m,order)
%    Euler_errors=accuracy(m,order,nxcuts)
%    Euler_errors=accuracy(m,order,nxcuts,nsims)
% 
%  Args:
% 
%     m (rise | dsge): scalar or vector of model objects.
% 
%     order ({solve_order} | integer): Approximation order for which to
%        compute the accuracy. Always less than or equal to solve_order
% 
%     nxcuts ({10} | integer): Number of quantiles per discretized
%       continuous shock
% 
%     nsims ({1000} | integer| 1x2 cell): Number of simulations to generate
%       the state of the endogenous variables. 
% 
% 
%  Returns:
%     :
% 
%     - **g** [cell | struct]: if the number of models is greater
%       that one, the output is a cell array of structures. Each structure
%       has the following fields:
%       - eqtn1, eqtn2,... : structures containing "mean" (mean of the euler
%       error across all simulations), "max" (maximum of the euler error
%       across all simulations) and "min" (minimum of the euler error across
%       all simulations)
% 
% 
%  Note:
% 
%     - The errors are NOT log10(abs(error)). The user can take the log10 or
%       any other transformation as desired.
% 
%     - The smaller the standard deviation of the shocks the higher the
%       accuracy even if the solution is not very accurate in the first
%       place.
% 
%     - The accuracy of the solution may critically depend on the tolerance
%       level used for solving the model in the first place. For instance,
%       if an iterative algorithm (e.g. mfi) is used, the solution might not
%       be 100% accurate even in a linear model
% 
%     - The errors are provided in absolute value for all equations, rather
%       than being normalized as is customarily done in the literature
%