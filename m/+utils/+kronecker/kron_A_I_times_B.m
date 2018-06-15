function C=kron_A_I_times_B(A,B,n)
% INTERNAL FUNCTION: computes C=kron(A,In)*B
%

[qn,m]=size(B);
[p,q]=size(A);
if qn/n~=q
    error('matrices sizes inconsistent')
end

C=reshape(reshape(B.',n*m,q)*A.',m,p*n).';
end
