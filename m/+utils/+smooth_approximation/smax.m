function a=smax(x,d)
% smooth approximation of max(x,0)
a=0.5*(x+utils.smooth_approximation.sabs(x,d));