function [B,BINT,R,RINT,STATS]=regress(this,this2,varargin)
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
% example: y=rise_time_series(1990,rand(100,1)); % random series
% X=y(-1)&y(-2)&y(-3); % columns of lags
% X=ones(X); % add a column of ones
% [B,BINT,R,RINT,STATS]=regress(y,X)
if isa(this,'rise_time_series') && isa(this2,'rise_time_series')
    [this,this2]=intersect(this,this2);
    y=double(this);
    if size(y,3)>1
        error([mfilename,':: this operation is only defined for databases with one page'])
    end
    if size(y,2)>1
        error([mfilename,':: first argument must have only one column of data'])
    end
    X=double(this2);
    if size(X,3)>1
        error([mfilename,':: this operation is only defined for databases with one page'])
    end
    [B,BINT,R,RINT,STATS]=regress(y,X,varargin{:});
else
    error([mfilename,':: both arguments should be ''rise_time_series'' objects'])
end
end
