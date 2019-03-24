%--- help for dsge/solve_alternatives ---
%
%  Attempts to find all solutions of a MSDSGE model
% 
%  ::
% 
%    [allobj,mean_square_stable] = solve_alternatives(obj)
%    [allobj,mean_square_stable] = solve_alternatives(obj,varargin)
% 
%  Args:
% 
%     obj (rise | dsge): RISE model object
%     solve_alternatives_nsim (numeric | {100}): Number of random starting
%       values generated.
% 
%  Returns:
%     :
% 
%     - **allobj** [rise|dsge]: vector of RISE objects with the solutions found
%