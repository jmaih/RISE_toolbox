function [y,x,n,T]=set_y_and_x(datay,datax,nlags,constant)
[n,smpl]=size(datay);
if constant
    datax=[datax;ones(1,smpl)];
end
nx=size(datax,1);
T=smpl-nlags;
y=datay(:,nlags+1:end);
x=nan(n*nlags+nx,T);
for ilag=1:nlags
    x((ilag-1)*n+1:n*ilag,:)=datay(:,(nlags+1:end)-ilag);
end
x(n*nlags+(1:nx),:)=datax(:,nlags+1:end);
