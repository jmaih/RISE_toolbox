%--- help for mdd/bridge ---
%
%  BRIDGE Computes the log marginal data density using the "bridge"
%  approximation by :cite:`MengWong96`. 
% 
%  ::
% 
%    log_mdd=BRIDGE(obj)
% 
%    log_mdd=BRIDGE(obj,fix_point)
% 
%    log_mdd=BRIDGE(obj,fix_point,opts)
% 
%  Args:
% 
%     - **obj** [mdd]: Marginal Data Density object
% 
%     - **fix_point** [{true}|false|empty]: if true, optimization is done
%       using a fix point algorithm. Else, an iterative strategy is used.
% 
%     - **opts** [empty|struct]: options see help for mdd.global_options
% 
%  Returns:
%     :
% 
%     - **log_mdd** [numeric]: log marginal data density
% 
%     - **extras** [struct]: history of log MDD and convergences at each
%       iteration 
% 
%  See also : mdd.laplace, mdd.mhm, mdd.is, mdd.ris, mdd.cj,
%    mdd.mueller, mdd.laplace_mcmc, mdd.swz, mdd.global_options
%