function A=ivech(v)
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


v=v(:);
nv=numel(v);
d=1-4*1*(-2*nv);
n=.5*(-1+sqrt(d));

A=zeros(n);
iter=0;
for ii=1:n
    A(ii:end,ii)=v(iter+(1:n-ii+1));
    iter=iter+n-ii+1;
end

A=A+transpose(tril(A,-1));