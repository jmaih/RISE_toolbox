function [h,p,jbstat,critval]=jbtest(this,varargin)
% H1 line
%
% ::
%
%
% Args:
%
% Returns:
%    :
%
% Note:
%
% Example:
%
%    See also:

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
this=this.data;
if size(this,3)>1
    error([mfilename,':: this operation is only defined for databases with one page'])
end
n=size(this,2);
h=nan(1,n);
p=nan(1,n);
jbstat=nan(1,n);
critval=nan(1,n);
for ii=1:n
    [h(ii),p(ii),jbstat(ii),critval(ii)]=jbtest(this(:,ii),varargin{:});
end
end
