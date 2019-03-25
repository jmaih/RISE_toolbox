%--- help for sym/laplace ---
%
% LAPLACE Laplace transform.
%    L = LAPLACE(F) is the Laplace transform of the sym F with default
%    independent variable t.  The default return is a function of s.
%    If F = F(s), then LAPLACE returns a function of z:  L = L(z).
%    By definition, L(s) = int(F(t)*exp(-s*t),t,0,inf).
% 
%    L = LAPLACE(F,z) makes L a function of z instead of the default s:
%    LAPLACE(F,z) <=> L(z) = int(F(t)*exp(-z*t),t,0,inf).
% 
%    L = LAPLACE(F,w,u) makes L a function of u instead of the
%    default s (integration with respect to w).
%    LAPLACE(F,w,u) <=> L(u) = int(F(w)*exp(-u*w),w,0,inf).
% 
%    Examples:
%     syms a s t w x F(t)
%     laplace(t^5)          returns   120/s^6
%     laplace(exp(a*s))     returns   -1/(a-z)
%     laplace(sin(w*x),t)   returns   w/(t^2+w^2)
%     laplace(cos(x*w),w,t) returns   t/(t^2+x^2)
%     laplace(x^(3/2),t)    returns   (3*pi^(1/2))/(4*t^(5/2))
%     laplace(diff(F(t)))   returns   s*laplace(F(t),t,s) - F(0)
% 
%    See also SYM/ILAPLACE, SYM/FOURIER, SYM/ZTRANS, SUBS.
%
%    Other functions named laplace
%
%       mdd/laplace
%