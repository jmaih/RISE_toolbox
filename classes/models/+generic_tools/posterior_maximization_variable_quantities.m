%  Recomputes the quantities that depend on the Hessian
% 
%  ::
% 
%    post_max=posterior_maximization_variable_quantities(post_max,H)
%    post_max=posterior_maximization_variable_quantities(post_max,H,flag)
% 
%  Args:
% 
%     - **post_max** [struct]: more specifically the content of
%       obj.estimation.posterior_maximization
%     - **a_func** [function_handle]: function that inflates x and Vx under
%       linear restrictions.
% 
%  Returns:
%     :
% 
%     - **post_max** [struct]: containing
% 
%       - **SD** [vector]: standard deviations of the parameters
%       - **log_marginal_data_density_laplace** [numeric]: laplace
%         approximation of the log marginal data density
%       - **vcov** [matrix]: variance-covariance of estimated parameters
% 
%