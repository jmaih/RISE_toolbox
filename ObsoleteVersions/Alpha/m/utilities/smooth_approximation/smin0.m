function a=smin0(x,d)
% smooth approximation of min(x,0)
a=0.5*(x-sabs(x,d));