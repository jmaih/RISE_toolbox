X=[];
const=1;
pex=0;
N=2;
p=4;
T=100;
parameters.N=N;
parameters.p=p;
parameters.Q=ones(const+(N*p)+pex,N)/1000;
parameters.A=[1,0;0.01,1];
parameters.Phi=ones(1,N)/100;
parameters.T=T;
shocks.epsB=randn(T,const+(N*p)+pex,N);
shocks.epsv=randn(T,N);
shocks.mu=randn(T,N);
parameters.B0=zeros(const+(N*p)+pex,N);
parameters.lambda0=0.01*ones(1,N);
[y,B,lambda,v]=svmodels(X,parameters,shocks,const);

% Ploting
figure
plot(y);
figure
plot([squeeze(B(:,:,1)),squeeze(B(:,:,2))]);
figure
plot(lambda);