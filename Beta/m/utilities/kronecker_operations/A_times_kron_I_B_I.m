function C=A_times_kron_I_B_I(A,B,q,r)
% computes C=A*kron(kron(Iq,B),Ir)

% clc,rb=10;cb=7;q=9;ra=3;r=4;ca=rb*q*r;
% B=rand(rb,cb);A=rand(ra,ca);
% max(max(abs(A_times_kron_I_B_I(A,B,q,r)-A*kron(kron(eye(q),B),eye(r)))))

if size(A,2)~=q*r*size(B,1)
    error('matrices sizes inconsistent')
end

C=A_times_kron_I_B(kron(B,eye(r)),A,q);

end
