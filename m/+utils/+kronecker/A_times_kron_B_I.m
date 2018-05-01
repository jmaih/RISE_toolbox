function C=A_times_kron_B_I(A,B,q,use_reshape)
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

% computes C=A*kron(B,Iq)

% clc,rb=10;cb=7;q=9;ra=3;ca=rb*q;
% B=rand(rb,cb);A=rand(ra,ca);
% max(max(abs(A_times_kron_B_I(A,B,q)-A*kron(B,eye(q)))))

if nargin<4
    use_reshape=true;
end

% sparse inputs
%---------------
if ~isa(A,'tsparse')
    A=sparse(A);
end
B=sparse(B);

[rb,cb]=size(B);
[ra,ca]=size(A);
if rb*q~=ca
    error('matrices sizes inconsistent')
end

if use_reshape
    C=reshape(A,ra*q,rb);
    C=C*B;
    C=reshape(C,ra,cb*q);
    % C=reshape(reshape(A,ra*q,rb)*B,ra,cb*q);
else
    C=A_times_kronBI();
end

% output is automatically sparse
%--------------------------------

    function C=A_times_kronBI()
        C=zeros(ra,cb*q);
        picka=0:q:(rb-1)*q;
        pickc=0:q:(cb-1)*q;
        for iq=q:-1:1
            C(:,pickc+iq)=A(:,picka+iq)*B;
        end
    end
end
