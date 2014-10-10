function C=A_times_kron_B_I(A,B,q)
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
% Description
% ------------
%
% Examples
% ---------
%
% See also: 

% computes C=A*kron(B,Iq)

% clc,rb=10;cb=7;q=9;ra=3;ca=rb*q;
% B=rand(rb,cb);A=rand(ra,ca);
% max(max(abs(A_times_kron_B_I(A,B,q)-A*kron(B,eye(q)))))

[rb,cb]=size(B);
[ra,ca]=size(A);
if rb*q~=ca
    error('matrices sizes inconsistent')
end

C=reshape(reshape(A,ra*q,rb)*B,ra,cb*q);
end
