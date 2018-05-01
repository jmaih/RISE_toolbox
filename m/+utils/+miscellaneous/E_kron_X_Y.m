function E=E_kron_X_Y(Mx,Ny,Oyx)
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

% problems:
%         1- X and Y should have the same # rows
%         2- How to get the COV(vec(Y),vec(X)) ?
%         3- can the double kronecker be avoided?

[nx,p]=size(Mx);

[ny,q]=size(Ny);

if nx~=ny
    error('matrices should have the same number of rows')
end

Ip=eye(p);

In=eye(n);

Kqn=utils.gridfuncs.commutation(q,n);

E=kron(kron(Ip,Kqn),In)*(vec(Oyx)+kron(vec(Mx),vec(Ny)));

E=reshape(E,n^2,p*q);

end

function x=vec(X)
x=X(:);
end