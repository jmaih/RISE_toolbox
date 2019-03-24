%--- help for generic/hessian ---
%
%  Computes the hessian of the model at a specific point
% 
%  ::
% 
%    [obj,H,issue] = hessian(obj)
% 
%    [obj,H,issue] = hessian(obj,x)
% 
%    [obj,H,issue] = hessian(obj,x,varargin)
% 
%  Args:
% 
%     obj (rise | dsge | rfvar | svar): model object
% 
%     x ([] | vector): vector at which one wants to compute the hessian
% 
%     varargin: additional optional inputs among which the most relevant
%       for estimation is:
% 
%       - **hessian_type** [{'fd'}|'opg']: The hessian is either computed by
%         finite differences (fd) or by outer-product-gradient (opg)
% 
%  Returns:
%     :
% 
%     - **obj** [rise|dsge|rfvar|svar]: model object containing the new hessian
%       in case the model was previously estimated.
% 
%     - **H** [d x d matrix]: hessian matrix
% 
%     - **issue** [char|'']: any issue encountered during the computation of
%       the hessian
% 
%
%    Other functions named hessian
%
%       sym/hessian
%