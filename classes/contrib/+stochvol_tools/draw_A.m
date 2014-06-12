function [Avec,A] = draw_A(y,Z,B,Lambdat,parameters,priors)   
%% Code to draw (A|B_t,Lambda_t)
% INPUT
% priors: priors for mean and variance of A, precisely:
% ch0:  mean prior for correlation parameters for stochastic volatilities
% (normal), e.g A_mean = zeros(N*(N-1)/2,1);
% ssc0:  variance prior for correlation parameters for stochastic
% volatilities (normal); e.g A_var = 10000*eye(N*(N-1)/2,N*(N-1)/2);
% data
% B: time varying B
% LambdaT: time-varying volatilities
% Parameters
% L: number of lags


%%
N=size(y,2);
L=parameters.lags;
evar = stochvol_tools.innov(y,Z,B',N,T,L);
ch0=priors.Amean;
ssc0==priors.Avar;
k = 0;
for si = 2:N,
    lhs = Lambdat(2:end,si).^.5;
    yhs=(evar(:,:)')./(lhs*ones(1,N));
    yr = yhs(:,si);
    xr = -yhs(:,1:si-1);
    j = k+1;
    k = si-1+k;
    ch(j:k,1) = stochvol_tools.bayesreg(ch0(j:k),ssc0(j:k,j:k),1,yr,xr);
end
Avec = ch;
A = stochvol_tools.chofac(N,Avec);
