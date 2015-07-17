function [h,p,jbstat,critval]=jbtest(this,varargin)
% -H: H=0 means null hypothesis ("the data are normally
% distributed") cannot be rejected at the 5% significance
% level. H=1 means null can be rejected at the 5% level.
% - P: p-value
% - JBSTAT: test statistic
% - CRITVAL: critical value for test
% JBTEST treats NaNs in X as missing values, and ignores them.
% Test Statistic:  JBSTAT = N*(SKEWNESS^2/6 + (KURTOSIS-3)^2/24),
%                   where N is the sample size and the kurtosis of
%                   the normal distribution is defined as 3.
%  H = JBTEST(X,ALPHA) performs the test at significance level
%  ALPHA
n=numel(this);
this=double(this);
if size(this,3)>1
    error([mfilename,':: this operation is only defined for databases with one page'])
end
h=nan(1,n);
p=nan(1,n);
jbstat=nan(1,n);
critval=nan(1,n);
for ii=1:n
    dd=this(:,ii);
    dd=dd(~isnan(dd));
    if isempty(dd)
        error([mfilename,':: no valid observations to compute the standard deviation in column ',int2str(ii)])
    end
    [h(ii),p(ii),jbstat(ii),critval(ii)]=jbtest(dd,varargin{:});
end
end
