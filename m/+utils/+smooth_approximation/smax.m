function a=smax(x,tol)
% smooth approximation of max(x,0)
if nargin<2
    tol=1e-5;
end
a=0.5*(x+utils.smooth_approximation.sabs(x,tol));
end