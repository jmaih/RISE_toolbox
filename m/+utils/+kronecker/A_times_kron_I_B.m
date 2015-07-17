function C=A_times_kron_I_B(A,B,q)
% H1 line
%
% Syntax
% -------
% ::
%
% Inputs
% -------
%
% Outputs
% --------
%
% More About
% ------------
%
% Examples
% ---------
%
% See also: 

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
iter_d=1:rb;
iter_c=1:cb;
for ii=1:q
    C(:,iter_c)=A(:,iter_d)*B;
    if ii<q
        iter_d=iter_d+rb;
        iter_c=iter_c+cb;
    end
end

% if at least one of the inputs is sparse, return sparse for we created a
% zeros (instead of sparse) matrix as output in order not to modify the
% structure of a sparse matrix with unknown nonzeros, which would be slow.
%--------------------------------------------------------------------------
if issparse(A) || issparse(B)
    C=sparse(C);
end

end