%--- help for mdd/ris ---
%
%  RIS Computes the log marginal data density using the "Reciprocal
%  Importance Sampling" approximation as in :cite:`FRUE2006`
% 
%  ::
% 
%    log_mdd=RIS(obj)
% 
%    log_mdd=RIS(obj,[],opts)
% 
%  Args:
% 
%     - **obj** [mdd]: Marginal Data Density object
% 
%     - **opts** [empty|struct]: options see help for mdd.global_options
% 
%  Returns:
% 
%     - **log_mdd** [numeric]: log marginal data density
% 
%  See also : mdd.laplace, mdd.mhm, mdd.is, mdd.cj,
%    mdd.mueller, mdd.bridge, mdd.laplace_mcmc, mdd.swz, mdd.global_options
%