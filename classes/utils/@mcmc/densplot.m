function hdl=densplot(obj,pname,chain_id,N)
% Plots marginal density of a given parameter
%
% ::
%
%    hdl=densplot(obj,pname)
%    hdl=densplot(obj,pname,chain_id)
%    hdl=densplot(obj,pname,chain_id,N)
%
% Args:
%
%    obj (mcmc object): mcmc object
%    pname (string): parameter name
%
%    chain_id (integer | {[]}): choice of chain for which to plot the
%      density. If empty, all chains are used.
%
%    N (integer | {250}): Number of point in the kernel density
%
% Returns:
%    :
%
%    - **hdl** [integer]: handle to the plot
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