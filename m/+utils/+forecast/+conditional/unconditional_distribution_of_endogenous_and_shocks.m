function [DYbar,OMG]=unconditional_distribution_of_endogenous_and_shocks(R,Yf)
OMG=R*R';
[ncv,ncp]=size(Yf);
DYbar=zeros(size(R,1),1);
DYbar(1:ncv*ncp)=Yf(:);
end
