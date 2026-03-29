%--- help for mdd/swz ---
%
%  SWZ Computes the log marginal data density using the :cite:`SimsWZ2008`
%  approximation
% 
%  ::
% 
%    log_mdd=SWZ(obj)
% 
%    log_mdd=SWZ(obj,swz_pvalue)
% 
%    log_mdd=SWZ(obj,swz_pvalue,opts)
% 
%  Args:
% 
%     - **obj** [mdd]: Marginal Data Density object
% 
%     -  **swz_pvalue** [{90}|empty|scalar]: scalar for the computation of
%        the lower bound. Must be greater than 0 and less than 100
% 
%     - **opts** [empty|struct]: options see help for mdd.global_options
% 
%  Returns:
%     :
% 
%     - **log_mdd** [numeric]: log marginal data density
% 
%     - **extras** [struct]: lower bound and corresponding quantile
% 
%  See also : mdd.laplace, mdd.mhm, mdd.is, mdd.ris, mdd.cj,
%    mdd.mueller, mdd.bridge, mdd.laplace_mcmc, mdd.global_options
%