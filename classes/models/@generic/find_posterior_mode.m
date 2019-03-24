%--- help for find_posterior_mode ---
%
%  INTERNAL FUNCTION: Finds the mode of the posterior
% 
%  Note:
%     - Though one can start MCMC from any point in theory, a good starting point with a good estimate of the covariance/Hessian is requirement to create a good proposal distribution for Bayesian analysis of DSGE models. This function finds the mode of the posterior for a candidate starting point of the Bayesian analysis.
%     - Though available for modularity of the RISE toolbox, one should use :func:`estimate <dsge.estimate>` for most use cases.
% 
%