%--- help for dsge/link_parameters ---
%
%  LINK_PARAMETERS Dynamically binds parameters together by creating links
%  between parameters for future evaluation. 
% 
%     m = LINK_PARAMETERS(m, expressions)
% 
%     This function is responsible for dynamically binding parameters
%     together by creating links between parameters for future evaluation. 
% 
%     - `m`: Scalar or vector of model objects. Each model may have multiple
%       parameterizations. 
%     - `expressions`: Character or cell array of strings representing the
%       expressions to bind. For example, expressions = 'alpha = beta + gamma'. 
% 
%     Returns:
%     - `m`: Updated model object with dynamically bound parameters.
% 
%     Example:
%          m = LINK_PARAMETERS(m, expressions)
% 
%     See also: dsge.view_linked_parameters, dsge.unlink_parameters
%