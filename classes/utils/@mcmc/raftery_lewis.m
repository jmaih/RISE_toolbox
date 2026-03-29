%--- help for mcmc/raftery_lewis ---
%
%  raftery_lewis computes :cite:`RafteryLewis1992` convergence diagnostics
% 
%  ::
% 
%     convergence=raftery_lewis(this)
%     convergence=raftery_lewis(this,Opts)
%     convergence=raftery_lewis(this,Opts,chain_id)
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