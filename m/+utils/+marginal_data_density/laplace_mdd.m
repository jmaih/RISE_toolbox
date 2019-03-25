%  Computes the log marginal data density using the laplace approximation
% 
%  ::
% 
%    log_mdd=laplace_mdd(log_post,Hinv)
% 
%  Args:
% 
%     - **log_post** [numeric]: log of posterior kernel evaluated at the mode
%     - **Hinv** [matrix]: inverse Hessian (negative of the second derivatives
%       of the log posterior kernel)
% 
%  Returns:
%     :
% 
%     - **log_mdd** [numeric]: log marginal data density
% 
%