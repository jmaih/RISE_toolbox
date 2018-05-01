function C=A_times_kron_I_B(A,B,q,loop_backward)
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

% computes C=A*kron(Iq,B)

% clc,rb=10;cb=7;q=9;ra=3;ca=rb*q;
% B=rand(rb,cb);A=rand(ra,ca);
% max(max(abs(A_times_kron_I_B(A,B,q)-A*kron(eye(q),B))))

if nargin<4
    loop_backward=false;
end

% sparse inputs
%---------------
A=sparse(A);
B=sparse(B);

[rb,cb]=size(B);
[ra,ca]=size(A);
if ca~=rb*q
    error('matrices sizes inconsistent')
end

C=zeros(ra,q*cb);

if loop_backward
    % no evidence that looping backward is faster
    %--------------------------------------------
    iter_d=(1:rb)+(q-1)*rb;
    iter_c=(1:cb)+(q-1)*cb;
    for ii=q:-1:1
        C(:,iter_c)=A(:,iter_d)*B;
        if ii>1
            iter_d=iter_d-rb;
            iter_c=iter_c-cb;
        end
    end
else
    iter_d=1:rb;
    iter_c=1:cb;
    for ii=1:q
        C(:,iter_c)=A(:,iter_d)*B;
        if ii<q
            iter_d=iter_d+rb;
            iter_c=iter_c+cb;
        end
    end
end

% sparse output
%--------------
C=sparse(C);

end