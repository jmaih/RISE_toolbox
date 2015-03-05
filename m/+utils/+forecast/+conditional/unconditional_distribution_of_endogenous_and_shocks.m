function [DYbar,OMG]=unconditional_distribution_of_endogenous_and_shocks(R,Yf)
% N.B: R contains information on both endogenous and exogenous

% The covariance repair step is expensive and probably not necessary
% OMG=utils.cov.nearest(R*R'); 
OMG=R*R'; 
[ncv,ncp]=size(Yf);
DYbar=zeros(size(R,1),1);
% endogenous come first
DYbar(1:ncv*ncp)=Yf(:);
% exogenous are forecasted to be zeros.
end
