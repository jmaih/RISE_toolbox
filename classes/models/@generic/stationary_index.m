%--- help for stationary_index ---
%
%  INTERNAL FUNCTION: Index for stationary variables
% 
%  ::
% 
%    ind=stationary_index(obj)
%    ind=stationary_index(obj,too_small)
% 
%  Args:
% 
%     obj (rise | dsge): model object
%     too_small (numeric | {1e-10}): cutoff criterion
% 
%  Returns:
%     :
% 
%     - **ind** [logical]: true if variable is stationary
% 
%