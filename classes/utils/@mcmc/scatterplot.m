%--- help for mcmc/scatterplot ---
%
%  Make a scatter plot from the mcmc chains
% 
%  ::
% 
%     hdl = scatterplot(mcobj,pname1,pname2);
%     hdl = scatterplot(mcobj,pname1,pname2,chain_id);
%     hdl = scatterplot(mcobj,pname1,pname2,chain_id,varargin);
% 
%  Args:
%     mcobj (mcmc object): mcmc object
%     pname1 (str): x-axis parameter of scatter plot
%     pname2 (str): y-axis parameter of scatter plot
%     chain_id (vector of int): id of the chain to use
%     varargin (varargin): options fed into `scatter <https://www.mathworks.com/help/matlab/ref/scatter.html>`_ function of matlab
% 
%  Returns:
%     :
% 
%     - **hdl** (handle object): handle to plot object
% 
%