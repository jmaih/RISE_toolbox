function K=commutation(m,n,fast)
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

% properties of the commutation matrix
%                i) Kmn'=Knm
%                ii) inv(Kmn)=Knm
%                iii) K1n=Kn1=In
%                iv) trace(Kmn)=1+d(m-1,n-1) where d(x,y) is the greatest
%                common divisor of x and y
%                v) det(Kmn)=(-1)^(0.25*m*(m-1)*n*(n-1))
%                vi) Kmn*vec(A)=vec(A')
%                vii) Kmn*(kron(A,B))*Kst =kron(B,A), where A is (n,s) and
%                B is (m,t). Equivalently, Kmn*(kron(A,B)) =kron(B,A)*Kts
%                viii) Let P be a matrix with m rows, Q a matrix with n
%                      rows, x an m-vector and y an n-vector and z
%                      arbitrary. We have
%                      kron(kron(z',P),y)=Kmn*kron(y*z',P)
%                      kron(kron(x,Q),z')=Kmn*kron(Q,x*z')
%                ix) trace(Kmn*kron(A',B)=trace(A'B)=vec(A')'*Kmn*vec(B)
%                where A and B are (m,n) matrices
if nargin<3
    fast=[];
    if nargin<2
        n=[];
    end
end
if isempty(fast)
    fast=true;
end
if isempty(n)
    n=m;
end

K=zeros(m*n);

if fast
    if fast==2
        nm=n*m;
        J=(1:n).';
        IND=nan(nm,2);
        offset=0;
    end
else
    proto=zeros(n,m);
end
for ii=1:m
    if fast==2
        IND(offset+J,:) = [(ii-1)*n+J,(J-1)*m+ii];
        offset=offset+n;
    else
        for jj=1:n
            if fast
                K((ii-1)*n+jj,(jj-1)*m+ii)=1;
            else
                proto_ij=proto;
                proto_ij(jj,ii)=1;
                K((ii-1)*n+(1:n),(jj-1)*m+(1:m))=proto_ij;
            end
        end
    end
end
if fast==2
    siz=nm*ones(1,2);
    IND = sub2ind(siz,IND(:,1),IND(:,2));
    K(IND)=1;
end

end