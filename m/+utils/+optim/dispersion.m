function d=dispersion(X,lb,ub)
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

ul=sqrt(eps)+ub-lb;
X=sort(X,2);
d=max(abs(X(:,1)-X(:,end))./ul);