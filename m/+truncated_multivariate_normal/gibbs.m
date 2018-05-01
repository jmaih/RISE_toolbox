function x=gibbs(mu,SIG,lb,ub,Nsim,burn)
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

if nargin<6
    burn=0;
    if nargin<5
        Nsim=1;
    end
end
npar=numel(mu);

x0=lb+rand(npar,1).*(ub-lb);
bad=isnan(x0);
x0(bad)=mu(bad); %
V=SIG\eye(npar);

x=nan(npar,Nsim);
for isim=1:Nsim+burn
    for k=1:npar
        j=setdiff(1:npar,k);
        Vjj=V(j,j)-V(j,k)*V(k,j)/V(k,k); % <--- Vjj=SIG(j,j)\eye(npar-1); % huge computational savings
        muk=mu(k)+SIG(k,j)*Vjj*(x0(j)-mu(j));
        sigk=sqrt(SIG(k,k)-SIG(k,j)*Vjj*SIG(j,k));
        x0(k)=muk+sigk*truncated_multivariate_normal.univariate_draw((lb(k)-muk)/sigk,(ub(k)-muk)/sigk);
    end
    if isim>burn
        x(:,isim-burn)=x0;
    end
end

end