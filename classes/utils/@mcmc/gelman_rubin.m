%--- help for gelman_rubin ---
%
%  gelman_rubin : computes...
% 
%  ::
% 
%     [p,R,W,B,Rp]=gelman_rubin(obj)
% 
%  Args:
% 
%     obj (mcmc object): mcmc object
% 
%     recursive (true|{false}):
% 
%  Returns:
%     :
% 
%     - **p** (struct): with parameter names as fields and for each
%     field a structure containing information on PSRF, Within and Between
%     variances and the variance. N.B: One of the fields name is
%     "multivariate_" and it represents the aggregated statistics.
% 
%     - **R** (Matrix): Potential scale reduction factor for each parameter
% 
%     - **W** (Matrix): Within variance
% 
%     - **B** (matrix): Between variance
% 
%     - **Rp** (vector): Multivariate Potential Scale Reduction Factor
% 
%  Warning:
% 
%     - This function requires multiple chains of MCMC samples. See
%       **nchain** option of samplers.
% 
%  References:
% 
%     - :cite:`gelman1992inference`
% 
%