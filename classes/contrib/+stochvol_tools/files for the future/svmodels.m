function [y,B,lambda,v]=svmodels(X,parameters,shocks,const)
%% Description
% Inputs: 
% X exogenous variables
% Parameters
% Q, Phi are variance of shocks for B and v
% A^{-1} covariance matrix
% Beta0 and lambda0 initial condition
% p lags
% pex exogenous lags
% N(0,1) Shocks
% epsB for B; epsv for v; mu for lambda
% Outputs:
% y simulated series
% B simulated tvp parameters
% lambda simulated sv (variance)
% mu simulated covariance matrix
% boolend
% const =1 for constant
%% Inputs
% parameters
Q=parameters.Q;
A=parameters.A;
invA=inv(A);
Phi=parameters.Phi;
p=parameters.p;
T=parameters.T;  % time lenght
N=parameters. N; % number of series 
pex=size(X,2);
% shocks
epsB=shocks.epsB;
epsv=shocks.epsv;
mu=shocks.mu;
% initial conditions
B0=parameters.B0; 
lambda0=parameters.lambda0; 

%% Model
y=zeros(T+p,N);
B=NaN(T+p,const+(N*p)+pex,N);
lnlambda=NaN(T+p,N);
v=NaN(T+p,N);
lnlambda(p,:)=log(lambda0);
B(p,:,:)=B0;
B(1:p-1,:,:)=zeros(p-1,const+(N*p)+pex,N);
y(p,:)=[const,reshape(y(p:-1:1,:),1,p*N),zeros(1,pex)]*B0;
for t=p+1:T+p
   lnlambda(t,:)=lnlambda(t-1,:)+sqrt(Phi).*mu(t-p,:);
   v(t,:)=invA*chol(diag(exp(lnlambda(t,:))))*epsv(t-p,:)';
   B(t,:,:)=squeeze(B(t-1,:,:))+sqrt(Q).*squeeze(epsB(t-p,:,:));
   y(t,:)=[const,reshape(y(t-1:-1:t-p,:),1,p*N)]*squeeze(B(t,:,:))+v(t,:);  % to add X
end
y=y(p+1:end,:);
B=B(p+1:end,:,:);
lambda=exp(lnlambda(p+1:end,:));
v=v(p+1:end,:);


