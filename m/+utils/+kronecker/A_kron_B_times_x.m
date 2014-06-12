function z=A_kron_B_times_x(A,B,x,xr,xc)
% computes the kron(A,B)*x
if nargin<4
    xc=size(A,2);
    if nargin<3
        xr=size(B,2);
    end
end
z=B*reshape(x,xr,xc)*transpose(A);
z=z(:);
end