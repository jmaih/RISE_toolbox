%--- help for dsge/resid ---
%
%  Compute the residuals from the steady state
% 
%  ::
% 
%    r=resid(obj)
%    r=resid(obj,trim)
% 
%  Args:
% 
%     - obj (rise | dsge): model object
%     - trim (true|{false}|1x1 cell|1x2 cell):
%       * if false, print all residuals as is and Do not print equations
%       * if true, set to zero residuals that are less than a tolerance
%       level whose default is 1e-13. Do not print equations
%       * If it is a cell array, then it must be either {tol,false} or
%       {tol,true}. Either element can be empty, the default is {1e-13,false}
% 
%        - the first element is the tolerance level below which all
%          residuals are set to 0
%        - the second element is either "true" or "false" and decides
%          whether to additionally print the equation
% 
%  Returns:
%     :
% 
%     - **r** [vector|matrix]: residuals
% 
%  Note:
% 
%     - if no output is requested, the residuals are printed on screen
% 
%     - if the model has not been solved, resid will call sstate to try to
%       solve for the steady state, with option imposed set to true in order
%       to avoid reoptimization.
% 
%     - The trim option is for display only. If the function is called with
%       an output argument, trim is not used.
%