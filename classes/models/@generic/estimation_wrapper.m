%--- help for estimation_wrapper ---
%
%  environment for the estimation of parameters
% 
%  ::
% 
%    [x1,f1,H,issue,viol,obj,funevals]=estimation_wrapper(obj,action,x0,lb,ub,funevals)
% 
%  Args:
% 
%     obj (rise | dsge | svar | rfvar): model object
% 
%     action ('draw' | {'eval'} | 'estimate'): intended action
% 
%     x0 ([] | d x 1 vector): paramter vector for start of the optimization
%       or for log-posterior-kernel evaluation
% 
%     lb (d x 1 vector): lower bound of the parameters to estimate
% 
%     ub (d x 1 vector): upper bound of the parameters to estimate
% 
%     funevals ([] | numeric): function evaluations
% 
%  Returns:
%     :
% 
%     - **x1** [[]|d x 1 vector]: final paramter vector
% 
%     - **f1** [numeric]: value of objective function evaluated at x1
% 
%     - **H** [d x d x 2 array]: Hessian, such that H(:,:,1) is the hessian
%       computed returned from the optimizer and H(:,:,2) is the hessian computed
%       numerically by RISE.
% 
%     - **issue** [char|cellstr]: list of problems encountered during the
%       process
% 
%     - **viol** [[]|vector]: restriction violations
% 
%     - **obj** [rise|dsge|svar|rfvar]: model object possibly modified
% 
%     - **funevals** [[]|numeric]: function evaluations
% 
%