%--- help for problem_reduction ---
%
%  INTERNAL FUNCTION: Prepares the solution for checking Mean Square
%  Stability of the system
% 
%  ::
% 
%    [T,Q,n,h]=problem_reduction(obj)
% 
%  Args:
% 
%     obj (dsge | rise): solved model object
% 
%  Returns:
%     :
% 
%     - **T** [1 x h cell]: Companion form for different regimes. Each cell
%       contains an n x n matrix
%     - **Q** [h x h matrix]: Transition matrix
%     - **n** [integer]: number of endogenous variable
%     - **h** [integer]: number of regimes
% 
%  See also:
%     - svar/problem_reduction
% 
%
%    Other functions named problem_reduction
%
%       dsge/problem_reduction    generic/problem_reduction
%