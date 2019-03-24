%--- help for dsge/accuracy ---
%
%  Solves dsge models
% 
%  ::
% 
%    Euler_errors=accuracy(m)
%    Euler_errors=accuracy(m,nxcuts)
%    Euler_errors=accuracy(m,nxcuts,nsims)
%    Euler_errors=accuracy(m,nxcuts,nsims,nburn)
%    Euler_errors=accuracy(m,nxcuts,nsims,nburn,reset_params)
%    [Euler_errors,outsims]=accuracy(...)
% 
%  Args:
% 
%     m (rise | dsge): scalar or vector of model objects.
% 
%     nxcuts ({10} | integer): Number of quantiles per discretized
%       continuous shock
% 
%     nsims ({1000} | integer| 1x2 cell): Number of simulations to generate
%       the state of the endogenous variables. Alternatively, if the input
%       is a cell, no simulations will be performed. In this case nsims is
%       of the form nsims={Y,regs}, where:
%       - Y : simulations of endogenous variables
%       - regs : corresponding simulations for regimes
% 
%     nburn ({100} | integer): Number of burned simulations
% 
%     reset_params ({true} | false): reset the parameters to their correct
%       values in evaluating the model equations.
% 
%  Returns:
%     :
% 
%     - **Euler_errors** [cell | struct]: if the number of models is greater
%       that one, the output is a cell array of structures. Each structure
%       has the following fields:
%       - EE [1 x nsols cell] : Euler equation errors for all the solutions
%       with each cell for one solution
%       - eqtn1, eqtn2,... : structures containing "mean" (mean of the euler
%       error across all simulations), "max" (maximum of the euler error
%       across all simulations) and "min" (minimum of the euler error across
%       all simulations)
% 
%     - **outsims** [1 x 2 x nsols cell]: Simulations for the endogenous
%       variables (first element) and the regimes (second element) for each
%       solution (3rd dimension)
% 
%  Note:
% 
%     - The errors are log10(abs(error))
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
%