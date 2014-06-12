function [T,eigvals]=companion_form(B,nlags,endo_nbr,nz)
% the VAR has the form y=A1*y{-1}+A2*y{-2}+...+Ap*y{-p}+C*z+u
dd=endo_nbr*(nlags-1);
T=[B(:,1:end-nz);
    eye(dd),zeros(dd,endo_nbr)];
if nargout>1
    eigvals=abs(eig(T));
end
end
