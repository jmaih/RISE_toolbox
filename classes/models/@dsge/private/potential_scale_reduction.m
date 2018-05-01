function [R,B,W]=potential_scale_reduction(x)
% H1 line
%
% ::
%
%
% Args:
%
% Returns:
%    :
%
% Note:
%
% Example:
%
%    See also:


[npar,Nsim]=size(x);
R=nan(npar,Nsim-1);
B=nan(npar,Nsim-1);
W=nan(npar,Nsim-1);

for ii=1:Nsim-1
    [R(:,ii),B(:,ii),W(:,ii)]=sequential_potential_scale_reduction(x(:,1:ii+1));
end

function [R,B,W]=sequential_potential_scale_reduction(x)
[m,n]=size(x);
if m<1
    error([mfilename,':: potential scale reduction factor requires at least 2 chains'])
end
xbarj=mean(x,2);
xbar=mean(xbarj);
s2j=1/(n-1)*sum(bsxfun(@minus,x,xbarj).^2,2);
B=n/(m-1)*sum((xbarj-xbar).^2);
W=mean(s2j);
V=(n-1)/n*W+1/n*B;
R=sqrt(V/W);