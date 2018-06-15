function d=dispersion(X,lb,ub)
% INTERNAL FUNCTION
%

ul=sqrt(eps)+ub-lb;
X=sort(X,2);
d=max(abs(X(:,1)-X(:,end))./ul);