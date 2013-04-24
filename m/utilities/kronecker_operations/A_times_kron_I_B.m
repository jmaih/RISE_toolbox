function C=A_times_kron_I_B(A,B,q)
% computes C=A*kron(Iq,B)

% clc,rb=10;cb=7;q=9;ra=3;ca=rb*q;
% B=rand(rb,cb);A=rand(ra,ca);
% max(max(abs(A_times_kron_I_B(A,B,q)-A*kron(eye(q),B))))

[rb,cb]=size(B);
[ra,ca]=size(A);
if ca~=rb*q
    error('matrices sizes inconsistent')
end

C=zeros(ra,q*cb);
for ii=1:q
    iter_d=(ii-1)*rb+1:ii*rb;
    iter_c=(ii-1)*cb+1:ii*cb;
    C(:,iter_c)=A(:,iter_d)*B;
end

end