function hdl=autocorrplot(obj,pname,chain_id,order)
% Make a plot of autocorrlation at different lags from the mcmc chains
%
% ::
%
%    hdl = autocorrplot(mcobj,pname);
%    hdl = autocorrplot(mcobj,pname,chain_id);
%    hdl = autocorrplot(mcobj,pname,chain_id,order);
%
% Args:
%    mcobj (mcmc object): mcmc object
%    pname (str): parameter to make the autocorrlation plot
%    chain_id (vector of int): id of the chain to use
%    order (int): maximum number of lags for the autocorrelation
%
% Returns:
%    :
%
%    - **hdl** (handle object): handle to plot object
%

if nargin<4

    order=[];

    if nargin<3

        chain_id=[];

    end

end

if isempty(order)

    order=40;

end

x=load_draws(obj,pname,chain_id);

nc=size(x,1);

autocorr=zeros(nc,order+1);

for io=1:order+1

    ii=io-1;

    for ic=1:nc

        xx=x(ic,ii+1:end).';

        yy=x(ic,1:end-ii).';

        autocorr(ic,io)=corr(xx,yy);

    end

end

start_at=2;

select=start_at:order;

if ic>1

    select=utils.plot.locate_ticks(order,5,start_at);

end

autocorr=autocorr(:,select);

hdl0=bar(select(:),autocorr.');

% hdl0=plot(autocorr(:,2:end).');

title(pname)

axis tight

if nargout

    hdl=hdl0;

end

end
