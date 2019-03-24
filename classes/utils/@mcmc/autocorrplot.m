%--- help for mcmc/autocorrplot ---
%
%  Plots autocorrelations of a given parameter
% 
%  ::
% 
%     hdl=autocorrplot(obj,pname)
%     hdl=autocorrplot(obj,pname,chain_id)
%     hdl=autocorrplot(obj,pname,chain_id,order)
% 
%  Args:
% 
%     obj (mcmc object): mcmc object
% 
%     pname (string): parameter name
% 
%     chain_id (integer | {[]}): choice of chain for which to plot the
%       autocorrelation. If empty, all chains are used.
% 
%     order (integer | {40}): maximum order of the autocorrelation
% 
%  Returns:
%     :
% 
%     - **hdl** [integer]: handle to the plot
% 
%