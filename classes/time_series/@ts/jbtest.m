function [h,p,jbstat,critval]=jbtest(this,varargin)
% JBTEST overloads Matlab's jbtest for ts objects. performs the Jarque-Bera
% goodness-of-fit test of composite normality 
%
% ::
%
%    [h,p,jbstat,critval]=jbtest(this)
%    [h,p,jbstat,critval]=jbtest(this,alpha)
%    [h,p,jbstat,critval]=jbtest(this,alpha,mctol)
%
% Args:
%
%    - **this** [ts]: time series object
%
%    - **alpha** [numeric|{0.05}]: significance level
%
%    - **mctol** [numeric|{[]}]: significance level when Monte Carlo is
%      used (instead of interpolation) in the computation of **p** (see
%      below)
%
% Output Args:
%
%    - **h** [0|1]: result of the test. H=0 indicates that the null
%      hypothesis ("the data are normally distributed") cannot be rejected
%      at the 5% significance level. H=1 indicates that the null hypothesis
%      can be rejected at the 5% level. 
%
%    - **p** : p-value computed using inverse interpolation into the
%      look-up table of critical values. Small values of p cast doubt on
%      the validity of the null hypothesis 
%
%    - **jbstat** : test statistic
%
%    - **critval** : critical value for the test
%
% See also : JBTEST

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
