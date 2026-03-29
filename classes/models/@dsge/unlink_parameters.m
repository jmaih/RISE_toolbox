%--- help for dsge/unlink_parameters ---
%
%  UNLINK_PARAMETERS remove the selected links 
% 
%     m = UNLINK_PARAMETERS(m, expressions)
% 
%     This function is responsible for unbinding dynamically bound parameters
% 
%     - `m`: Scalar or vector of model objects. 
%     - `targets`: Character or cell array of strings representing the
%       parameters to unlink. e.g. targets = {'alpha','beta','gamma'} will remove
%       the definitions of the specified parameters. It could also be '*',
%       in which case all links are deleted.
% 
%     Returns:
%     - `m`: Updated model object with the remaining dynamically bound parameters.
% 
%     Example:
%          m = UNLINK_PARAMETERS(m, {'alpha','beta','gamma'})
% 
%     See also: dsge.view_linked_parameters, dsge.link_parameters
%