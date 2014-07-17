function a=smin(x,d)
% smooth approximation of min(x,0)
a=0.5*(x-utils.smooth_approximation.sabs(x,d));