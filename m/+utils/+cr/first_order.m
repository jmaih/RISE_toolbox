%  INTERNAL FUNCTION: first-order multivariate chain rule
% 
%  ::
% 
%    res=first_order(dv,vz)
% 
%  Args:
% 
%     - **dv** [nd x nv matrix]: jacobian of function with respect to the
%       locations of its arguments
%     - **vz** [nv x nz matrix]: jacobian of the locations with respect to the
%       variables to differentiate
% 
%  Returns:
%     :
% 
%     - **res** [nd x nz]: output matrix
% 
%