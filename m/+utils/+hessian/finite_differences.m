%  Computes the hessian by finite differences
% 
%  ::
% 
%     H = finite_differences(Objective,params)
%     H = finite_differences(Objective,params,varargin)
% 
%  Args:
% 
%     - **Objective** [char|function handle]: function to differentiate
%     - **params** [vector]: point at which the differentiation is taken
%     - **varargin** : optional/further arguments of the objective function
% 
%  Returns:
%     :
% 
%     - **H** [matrix]: Hessian matrix
% 
%  See also:
% 
%     - utils.hessian.outer_product
% 
%