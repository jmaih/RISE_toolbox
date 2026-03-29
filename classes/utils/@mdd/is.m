%--- help for mdd/is ---
%
%  IS Computes the log marginal data density using the "importance sampling"
%  approximation
% 
%  ::
% 
%    log_mdd=IS(obj)
% 
%    log_mdd=IS(obj,[],opts)
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
%  See also : mdd.laplace, mdd.mhm, mdd.ris, mdd.cj,
%    mdd.mueller, mdd.bridge, mdd.laplace_mcmc, mdd.swz, mdd.global_options
%
%    Other uses of is
%
%       matlab.ui.internal.controller.FigureController/is
%