function C=A_times_k_kron_B(A,B,k)
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


% computes C=A*(BxBx...xB)

% clc,rb=10;k=3;ra=5;ca=rb^k; cb=30;
% B=rand(rb,cb);A=rand(ra,ca);
% max(max(abs(A_times_k_kron_B(A,B,k)-A*kron(B,kron(B,B)))))

[rb,cb]=size(B);
[~,ca]=size(A);
if ca~=rb^k
    error('matrices sizes inconsistent')
end

C=utils.kronecker.A_times_kron_B_I(A,B,rb^(k-1));
for ii=3:k
    C=utils.kronecker.A_times_kron_I_B_I(C,B,cb^(ii-2)*rb^(k-ii),rb);
end
C=utils.kronecker.A_times_kron_I_B(C,B,cb^(k-1));

% C=utils.kronecker.A_times_kron_I_B_I(A,B,q,r) % cb=q*r*rc;
% computes C=A*kron(kron(Iq,B),Ir)
end