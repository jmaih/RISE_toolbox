%  INTERNAL FUNCTION: memoizer for some dirichlet routines
% 
%  ::
% 
%    d=dirichlet_shortcuts()
%    d=dirichlet_shortcuts(a,location)
%    d=dirichlet_shortcuts(a,location,lpdfn,rndfn)
%    d=dirichlet_shortcuts(a,location,lpdfn,rndfn,sum_aij)
% 
%  Args:
% 
%     - **a** [k x 1 vector]: of hyperparameters for the distribution
%     - **location** [k-1 x 1]: location of the parameters to estimate. Note
%       that this vector has k-1 elements instead of k
%     - **lpdfn** [function_handle]: function handle for the log pdf of the
%       dirichlet distribution
%     - **rndfn** [function_handle]: function handle for random draws of the
%       dirichlet distribution
%     - **sum_aij** [numeric]: sum of the weights including the diagonal term
% 
%  Returns:
%     :
% 
%     - **d** [struct]: with the following elements
% 
%       - **lpdfn** [function_handle]: takes as input an k-1 x 1 vector but
%         returns the log density for an k x 1 vector
%       - **rndfn** [function_handle]: returns and k-1 x n matrix of random
%         draws of the dirichlet distribution
%       - **location** [vector]: location of the parameters of interest in the
%         vector of estimated parameters
%       - **sum_aij** [numeric]: sum of the weights including the diagonal term
% 
%