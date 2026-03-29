%--- help for mcmc/geweke ---
%
%  geweke computes Geweke convergence diagnostics
% 
%  ::
% 
%     convergence=geweke(this)
%     convergence=geweke(this,Opts)
%     convergence=geweke(this,Opts,chain_id)
% 
%  Args:
% 
%     this (mcmc object): mcmc object
% 
%     pname (string): parameter name
% 
%     chain_id (integer | {1}): choice of chain for which to do the
%     convergence analysis
% 
%  Returns:
%     :
% 
%     - **convergence** [struct]: convergence analysis for each parameter
%