%--- help for mdd/mueller ---
%
%  MUELLER Computes the log marginal data density using the Ulrich Mueller's
%  approximation
% 
%  ::
% 
%    log_mdd=MUELLER(obj)
% 
%    log_mdd=MUELLER(obj,[],opts)
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
%  See also : mdd.laplace, mdd.mhm, mdd.is, mdd.ris, mdd.cj,
%    mdd.mueller, mdd.bridge, mdd.laplace_mcmc, mdd.swz, mdd.global_options
%