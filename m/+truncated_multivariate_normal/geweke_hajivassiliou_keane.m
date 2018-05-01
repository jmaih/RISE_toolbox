function x=geweke_hajivassiliou_keane(mu,CS,lb,ub,Nsim)
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
x=nan(npar,Nsim);
lb1=lb-mu;
ub1=ub-mu;
for k=1:npar
    vx=CS(k,1:k-1)*x(1:k-1,:);
    tmp_low=(lb1(k)*ones(1,Nsim)-vx)/CS(k,k);
    tmp_high=(ub1(k)*ones(1,Nsim)-vx)/CS(k,k);
    
    x(k,:)=truncated_multivariate_normal.univariate_draw(tmp_low,tmp_high);
end
x=mu(:,ones(1,Nsim))+CS*x;

end