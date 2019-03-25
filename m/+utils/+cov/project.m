%  INTERNAL FUNCTION: Projection of covariance matrix such the eigenvalues
%  are sandwiched. 
% 
%  ::
% 
%     vcov = project(vcov0,e_min,e_max);
% 
%  Args:
% 
%     - **vcov0** [matrix]: initial covariance matrix
%     - **e_min** [[]|{sqrt(eps)}]: scalar such that the minimum eigenvalue of vcov
%       is greater than or equal to "e_min".
%     - **e_max** [[]|{1/e_min}]: scalar such that maximum eigenvalue of vcov
%       is less than or equal to "e_max"
% 
%  Returns:
%     :
% 
%     - **vcov** [matrix]: updated covariance matrix
% 
%