function z=A_kron_B_times_x(A,B,x,cols_b,cols_a)
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

% computes the kron(A,B)*x
if nargin<4
    cols_a=size(A,2);
    if nargin<3
        cols_b=size(B,2);
    end
end
z=B*reshape(x,cols_b,cols_a)*transpose(A);
z=z(:);
end