function C=A_times_k_kron_B(A,B,k)

% computes C=A*(BxBx...xB)

% clc,rb=10;k=3;ra=5;ca=rb^k;
% B=rand(rb);A=rand(ra,ca);
% max(max(abs(A_times_k_kron_B(A,B,k)-A*kron(B,kron(B,B)))))

[~,cb]=size(B);
[~,ca]=size(A);
if ca~=cb^k
    error('mbtrices sizes inconsistent')
end

nk1=cb^(k-1);
nk2=cb^(k-2);
C=A_times_kron_B_I(A,B,nk1);
for ii=2:k-1
    C=A_times_kron_I_B_I(B,C,nk2,cb);
end
C=A_times_kron_I_B(B,C,nk1);

end