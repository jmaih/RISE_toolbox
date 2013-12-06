function a=sabs(x,d)
% smooth approximation of abs
a=x.^2./sqrt(x.^2+d.^2);