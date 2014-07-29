function a=sabs(x,tol)
if nargin<2
    tol=1e-5;
end
% smooth approximation of abs
a=x.^2./sqrt(x.^2+tol.^2);
end