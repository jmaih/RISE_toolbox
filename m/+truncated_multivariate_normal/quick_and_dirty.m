function x=quick_and_dirty(mu,CS,lb,ub,Nsim)
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

if nargin<5
    Nsim=1;
end
npar=numel(mu);

x=mu(:,ones(1,Nsim))+CS*randn(npar,Nsim);
lb=lb(:,ones(1,Nsim));
ub=ub(:,ones(1,Nsim));
xlow=x<lb;
x(xlow)=lb(xlow);
xhigh=x>ub;
x(xhigh)=ub(xhigh);

end