%--- help for mcmc/meanplot ---
%
%  Make a plot of cumulative means from the mcmc chain
% 
%  ::
% 
%     hdl = meanplot(mcobj,pname);
%     hdl = meanplot(mcobj,pname,chain_id);
% 
%  Args:
% 
%     mcobj (mcmc object): mcmc object
% 
%     pname (str): parameter to make the mean plot
% 
%     chain_id (vector of int): id of the chain to use
% 
%  Returns:
%     :
% 
%     - **hdl** (handle object): handle to plot object
% 
%