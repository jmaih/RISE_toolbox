%--- help for mcmc/psrf_plot ---
%
%  Make a plot of the "posterior scale reduction factor," i.e., Gelman-Rubin
%  diagnotics from the chains 
% 
%  ::
% 
%     hdl = psrf_plot(obj, pname)
% 
%  Args:
% 
%     obj (mcmc object): mcmc object
% 
%     pname (char): parameter name. N.B: One of the parameter names is
%     "multivariate_" and it represents the aggregated statistics.
% 
%     start (numeric|{1}|function handle): iteration at which to start the
%        plot of the PSRF. If a function handle is used then it should take
%        as input the total number of observations and return the point at
%        which to start. e.g. @(x)round(0.5*x)
% 
% 
%  Returns:
%     :
% 
%     - **hdl** (handle object): handle to plot object
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