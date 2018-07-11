function hdl=densplot(obj,pname,chain_id,N)
% Make a kernel density plot from the mcmc chains
%
% ::
%
%    hdl = densplot(mcobj,pname);
%    hdl = densplot(mcobj,pname,chain_id);
%    hdl = densplot(mcobj,pname,chain_id,N);
%
% Args:
%    mcobj (mcmc object): mcmc object
%    pname (str): parameter to make the density plot
%    chain_id (vector of int): id of the chain to use
%    N (int): number of sampling points (also determines the bandwith of the kernel)
%
% Returns:
%    :
%
%    - **hdl** (handle object): handle to plot object
%

if nargin<4

    N=[];

    if nargin<3

        chain_id=[];

    end

end

if isempty(N)

    N=250;

end

x=load_draws(obj,pname,chain_id);

x=vec(x.');

[f_kdens,x_kdens]=distributions.kernel_density(x,[],[],'normal',N);

hdl0=plot(x_kdens,f_kdens);

title(pname)

axis tight

if nargout

    hdl=hdl0;

end

end