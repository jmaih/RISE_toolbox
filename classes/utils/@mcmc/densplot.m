%--- help for mcmc/densplot ---
%
%  Plots marginal density of a given parameter
% 
%  ::
% 
%     hdl=densplot(obj,pname)
%     hdl=densplot(obj,pname,chain_id)
%     hdl=densplot(obj,pname,chain_id,N)
% 
%  Args:
% 
%     obj (mcmc object): mcmc object
% 
%     pname (string): parameter name
% 
%     chain_id (integer | {[]}): choice of chain for which to plot the
%       density. If empty, all chains are used.
% 
%     N (integer | {250}): Number of point in the kernel density
% 
%  Returns:
%     :
% 
%     - **hdl** [integer]: handle to the plot
% 
%