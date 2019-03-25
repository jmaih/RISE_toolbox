%  INTERNAL FUNCTION: generate candidates for optimization
% 
%  ::
% 
%    [x,f,viol,funevals]=generate_candidates(objective,lb,ub,n,...
%        restrictions,penalty,varargin)
% 
%  Args:
%     - **objective** [function_handle]: objective to minimize
%     - **lb** [vector]: lower bound of the search space
%     - **ub** [vector]: upper bound of the search space
%     - **n** [integer]: number of candidates to generate
%     - **max_trials** [integer]: number of trials after which the procedure
%       crashes
%     - **restrictions** [empty|function_handle]: function evaluating the
%       violations
%     - **opt** [struct]: structure with fields
% 
%       - **restrictions_in_objective** [true|false]:
%       - **returns_retcode** [true|false]:
%       - **restrictions_same_weights** [true|{false}]:
%       - **allow_restrictions_violations** [true|{false}]:
% 
%     - **penalty** [numeric]: value functions exceeding this value in
%       absolute value are assigned this value
%     - **varargin** []: additional input arguments for the objective function
% 
%  Returns:
%     :
% 
%     - **x** [d x n matrix]: parameter vectors
%     - **f** [row vector]: value function for each parameter vector
%     - **viol** [row vector]: violation penalty for each parameter vector
%     - **funevals** [integer]: number of function evaluations
% 
%