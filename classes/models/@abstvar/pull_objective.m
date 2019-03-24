%--- help for abstvar/pull_objective ---
%
%  Pulls the objective function to optimize
% 
%  ::
% 
%    [ff,lb,ub]=pull_objective(self)
%    [ff,lb,ub,x0]=pull_objective(self)
%    [ff,lb,ub,x0,vcov]=pull_objective(self)
%    [ff,lb,ub,x0,vcov,self]=pull_objective(self)
% 
%  Args:
% 
%     self (svar | rfvar): initial model object
% 
%  Returns:
%     :
% 
%     - **ff** [function handle]: for computing "minus log posterior kernel"
%     - **lb** [d x 1 vector]: lower bound of the parameters to optimize
%     - **ub** [d x 1 vector]: upper bound of the parameters to optimize
%     - **x0** [d x 1 vector]: posterior mode if available
%     - **vcov** [d x d matrix]: covariance matrix at the posterior mode if
%       available.
%     - **self** [rise|dsge|svar|rfvar]: updated model object
% 
%  Note:
% 
%     - The function can be used for :
% 
%        - optimization,
%        - gradient computation,
%        - hessian computation,
%        - posterior simulation
% 
%     - The updated object should be used for doing various exercises (irfs,
%       simulations, etc.) if the posterior mode is not computed.
% 
%     - Using this function is potentially costly, one could alternatively
%       simply use log_posterior_kernel. However, if there are restrictions, they
%       will not be enforced. Nevertheless it is an interesting proposition that
%       should be investigated further.
% 
%
%    Other functions named pull_objective
%
%       dsge/pull_objective    generic/pull_objective
%