%  Computes the hessian using the outer product of the gradient 
% 
%  ::
% 
%     H = outer_product(Objective,params);
%     H = outer_product(Objective,params,varargin);
% 
%  Args:
% 
%     - **Objective** [char|function handle]: function to differentiate
% 
%     - **params** [vector]: point at which the differentiation is taken
% 
%     - **varargin** : optional/further arguments of the objective function
% 
%  Returns:
%     :
% 
%     - **H** [matrix]: Hessian matrix
% 
%     - **H2** [matrix]: Hessian matrix scaled by the number of elements in
%        the underlying gradient
% 
%     - **grad** [vector]: gradient vector
% 
%  Note:
%     Objective is assumed to have at least two output arguments the
%     first one will not be used. The second one is the different increments of
%     the likelihood 
% 
%  See also:
% 
%     - utils.hessian.finite_differences
% 
%