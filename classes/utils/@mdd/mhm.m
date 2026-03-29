%--- help for mdd/mhm ---
%
%  mhm Computes the log marginal data density using the "modified harmonic
%  mean" approximation by :cite:`Geweke1999`.
% 
%  ::
% 
%    log_mdd=mhm(obj)
% 
%  Args:
% 
%     - **obj** [mdd]: Marginal Data Density object
% 
%     -  **mhm_tau** [{(.1:.1:.9)}|vector|scalar|empty]: truncation
%        probabilities 
% 
%     - **opts** [empty|struct]: options see help for mdd.global_options
% 
%  Returns:
% 
%     - **log_mdd** [numeric]: log marginal data density
% 
%  See also : mdd.laplace, mdd.is, mdd.ris, mdd.cj,
%    mdd.mueller, mdd.bridge, mdd.laplace_mcmc, mdd.swz, mdd.global_options
%