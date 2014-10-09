function [B,BINT,R,RINT,STATS]=ar(this,lag,const)
% - B: vector of regression coefficients in the linear model Y = X*B.
% - BINT: of 95% confidence intervals for B
% - R: vector of residuals
% - RINT: matrix of intervals that can be used to diagnose
% outliers.  If RINT(i,:) does not contain zero, then the i-th
% residual is larger than would be expected, at the 5%
% significance level.  This is evidence that the I-th
% observation is an  outlier.
% - STATS: vector containing, in the following order, the R-square
% statistic, the F statistic and p value for the full model,
% and an estimate of the error variance.
if nargin<3
    const=true;
    if nargin<2
        lag=1;
    end
end
Z=double(this);
if size(Z,3)>1
    error([mfilename,':: this operation is only defined for databases with one page'])
end
if size(Z,2)>1
    error([mfilename,':: first argument must have only one column of data'])
end
y=Z(lag+1:end);
smpl=size(y,1);
X=nan(smpl,lag);
for ii=1:lag
    X(:,ii)=Z((lag+1:end)-ii);
end
if const
    X=[X,ones(smpl,1)];
end
[B,BINT,R,RINT,STATS]=regress(y,X);
end
