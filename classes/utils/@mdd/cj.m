%--- help for mdd/cj ---
%
%  CJ Computes the log marginal data density using the
%  :cite:`ChibJeliazkov2001` approximation  
% 
%  ::
% 
%    log_mdd=CJ(obj)
% 
%    log_mdd=CJ(obj,[],opts)
% 
%  Args:
% 
%     - **obj** [mdd]: Marginal Data Density object
% 
%     - **opts** [empty|struct]: options see help for mdd.global_options
% 
%  Returns:
%     :
% 
%     - **log_mdd** [numeric]: log marginal data density
% 
%  See also : mdd.laplace, mdd.mhm, mdd.is, mdd.ris,
%    mdd.mueller, mdd.bridge, mdd.laplace_mcmc, mdd.swz, mdd.global_options
%