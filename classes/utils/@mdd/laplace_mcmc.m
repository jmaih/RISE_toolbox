%--- help for mdd/laplace_mcmc ---
%
%  laplace_mcmc Computes the log marginal data density using the laplace
%  approximation and a covariance matrix computed from the draws of
%  posterior sampling
% 
%  ::
% 
%    log_mdd=laplace_mcmc(obj)
% 
%  Args:
% 
%     - **obj** [mdd]: Marginal Data Density object
% 
%  Returns:
%     :
% 
%     - **log_mdd** [numeric]: log marginal data density
% 
%  See also : mdd.laplace, mdd.mhm, mdd.is, mdd.ris, mdd.cj,
%    mdd.mueller, mdd.bridge, mdd.swz
%