function hdl=scatterplot(obj,pname1,pname2,chain_id,varargin)
% Make a scatter plot from the mcmc chains
%
% ::
%
%    hdl = scatterplot(mcobj,pname1,pname2);
%    hdl = scatterplot(mcobj,pname1,pname2,chain_id);
%    hdl = scatterplot(mcobj,pname1,pname2,chain_id,varargin);
%
% Args:
%    mcobj (mcmc object): mcmc object
%    pname1 (str): x-axis parameter of scatter plot
%    pname2 (str): y-axis parameter of scatter plot
%    chain_id (vector of int): id of the chain to use
%    varargin (varargin): options fed into `scatter <https://www.mathworks.com/help/matlab/ref/scatter.html>`_ function of matlab
%
% Returns:
%    :
%
%    - **hdl** (handle object): handle to plot object
%


if nargin<4

    chain_id=[];

end

x1=load_draws(obj,pname1,chain_id);

x1=vec(x1.');

x2=load_draws(obj,pname2,chain_id);

x2=vec(x2.');

hdl0=scatter(x1,x2,varargin{:});

xlabel(pname1)

ylabel(pname2)

axis tight

if nargout

    hdl=hdl0;

end

end