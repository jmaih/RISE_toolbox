function [DYbar,OMG]=UnconditionalDistributionOfEndogenousAndShocks(R,Yf)
OMG=R*R';
[ncv,ncp]=size(Yf);
DYbar=zeros(size(R,1),1);
DYbar(1:ncv*ncp)=Yf(:);
end
