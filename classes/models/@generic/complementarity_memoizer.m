%--- help for complementarity_memoizer ---
%
%  Memoize linear and nonlinear restrictions
% 
%  ::
% 
%    [sep_cf,cf]=complementarity_memoizer(obj)
% 
%  Args:
% 
%     obj (rise | dsge | svar | rfvar): model object
% 
%  Returns:
%     :
% 
%     - **sep_cf** [function handle]: function that takes as input a vector of
%       variable values and returns a vector of values that are expected to be
%       greater than or equal to 0.
% 
%     - **cf** [function handle]: function that takes as input a vector of
%       variable values and returns a true if all the restrictions are satisfied
%       and returns false otherwise
% 
%