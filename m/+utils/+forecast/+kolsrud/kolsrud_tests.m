clc
N=10000;
T=20;
G=35;
Y=randn(N,T,G);
c=utils.forecast.kolsrud.chebyshev_distance(Y);
mvcb=utils.forecast.kolsrud.multivariate_chebyshev_box(Y,0.68);
size(mvcb)
