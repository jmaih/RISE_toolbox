function C=A_times_kron_I_B_I(A,B,q,r)
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

% computes C=A*kron(kron(Iq,B),Ir)

% clc,rb=10;cb=7;q=9;ra=3;r=4;ca=rb*q*r;
% B=rand(rb,cb);A=rand(ra,ca);
% max(max(abs(A_times_kron_I_B_I(A,B,q,r)-A*kron(kron(eye(q),B),eye(r)))))

[rb,cb]=size(B);
ca=q*r*rb;
if size(A,2)~=ca
    error('matrices sizes inconsistent')
end

C=zeros(size(A,1),q*r*cb);
rrb=r*rb;
rcb=r*cb;
for icol=1:q
    iterc=(icol-1)*rcb+1:icol*rcb;
    itera=(icol-1)*rrb+1:icol*rrb;
    C(:,iterc)=utils.kronecker.A_times_kron_B_I(A(:,itera),B,r);%A(:,itera)*kron(B,eye(r))
end

end
