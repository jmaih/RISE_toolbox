function hdl=meanplot(obj,pname,chain_id)
% Make a plot of cumulative means from the mcmc chain
%
% ::
%
%    hdl = meanplot(mcobj,pname);
%    hdl = meanplot(mcobj,pname,chain_id);
%
% Args:
%
%    mcobj (mcmc object): mcmc object
%
%    pname (str): parameter to make the mean plot
%
%    chain_id (vector of int): id of the chain to use
%
% Returns:
%    :
%
%    - **hdl** (handle object): handle to plot object
%

if nargin<3

    chain_id=[];

end

x=load_draws(obj,pname,chain_id);

m=x;

nc=size(x,1);

V=zeros(nc,nc,obj.npop);

for ii=2:obj.npop

    [m(:,ii),V(:,:,ii)]=utils.moments.recursive(m(:,ii-1),V(:,:,ii-1),x(:,ii),ii);

end

t=(1:obj.npop)+obj.start-1;

hdl0=plot(t.',m.');

title(pname)

axis tight

if nargout

    hdl=hdl0;

end

end