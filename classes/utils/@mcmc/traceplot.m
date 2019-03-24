%--- help for mcmc/traceplot ---
%
%  Make a trace plot from the mcmc chains
% 
%  ::
% 
%     hdl = traceplot(mcobj,pname);
%     hdl = traceplot(mcobj,pname,chain_id);
%     hdl = traceplot(mcobj,pname,chain_id,ma_window);
% 
%  Args:
% 
%     mcobj (mcmc object): mcmc object
% 
%     pname (str): parameter to make the trace plot
% 
%     chain_id (vector of int): id of the chain to use
% 
%     ma_window (int): window size if using moving average smoothing.
% 
%  Returns:
%     :
% 
%     - **hdl** (handle object): handle to plot object
% 
%