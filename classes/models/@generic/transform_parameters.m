%--- help for transform_parameters ---
%
%  transform the parameters under linear restrictions and or dirichlet priors
% 
%  ::
% 
%    [obj,x0,lb,ub,vcov]=transform_parameters(obj,x0,lb,ub)
%    [obj,x0,lb,ub,vcov]=transform_parameters(obj,x0,lb,ub,vcov)
% 
%  Args:
% 
%     obj (rise | dsge | rfvar | svar): model object
% 
%     x0 (empty | n x 1 vector): initial conditions of estimation in a space that
%       the user understands
% 
%     lb (n x 1 vector): lower bound of the search space
% 
%     ub (n x 1 vector): upper bound of the search space
% 
%     vcov (empty | n x n matrix): variance covariance of the parameters
% 
%  Returns:
%     :
% 
%     - **obj** [rise|dsge|rfvar|svar]: model object with all restrictions set
%       up
% 
%     - **x0** [empty|m x 1 vector]: transform initial conditions with m<=n
% 
%     - **lb** [m x 1 vector]: transformed lower bound of the search space with
%       m<=n
% 
%     - **ub** [m x 1 vector]: transformed upper bound of the search space with
%       m<=n
% 
%     - **vcov** [empty|m x m matrix]: transformed variance covariance of the
%       parameters with m<=n
% 
%  Note:
% 
%     - Some checks have to be made after transformation in order to insure
%       that no element in the transformed lower bound exceeds its upper bound
%       counterpart.
% 
%  See also:
%     - untransform_parameters
% 
%