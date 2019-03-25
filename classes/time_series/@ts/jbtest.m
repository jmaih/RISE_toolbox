% JBTEST Jarque-Bera hypothesis test of composite normality.
%    H = JBTEST(X) performs the Jarque-Bera goodness-of-fit test of composite
%    normality, i.e., that the data in the vector X came from an unspecified
%    normal distribution, and returns the result of the test in H. H=0
%    indicates that the null hypothesis ("the data are normally distributed")
%    cannot be rejected at the 5% significance level. H=1 indicates that the
%    null hypothesis can be rejected at the 5% level.
% 
%    JBTEST treats NaNs in X as missing values, and ignores them.
% 
%    The Jarque-Bera test is a 2-sided goodness-of-fit test suitable for
%    situations where a fully-specified null distribution is not known, and its
%    parameters must be estimated.  For large sample sizes, the test statistic
%    has a chi-square distribution with two degrees of freedom. Critical
%    values, computed using Monte-Carlo simulation, have been tabulated for
%    sample sizes N <= 2000 and significance levels 0.001 <= ALPHA <= 0.50.
%    JBTEST computes a critical value for a given test by interpolating into
%    that table, using the analytic approximation to extrapolate for larger
%    sample sizes.
% 
%    The Jarque-Bera hypotheses and test statistic are:
% 
%              Null Hypothesis:  X is normally distributed with unspecified
%                                mean and standard deviation.
%       Alternative Hypothesis:  X is not normally distributed.  The test is
%                                specifically designed for alternatives in the
%                                Pearson family of distributions.
%               Test Statistic:  JBSTAT = N*(SKEWNESS^2/6 + (KURTOSIS-3)^2/24),
%                                where N is the sample size and the kurtosis of
%                                the normal distribution is defined as 3.
% 
%    H = JBTEST(X,ALPHA) performs the test at significance level ALPHA.  ALPHA
%    is a scalar in the range 0.001 <= ALPHA <= 0.50.  To perform the test at
%    significance levels outside that range, use the MCTOL input argument.
% 
%    [H,P] = JBTEST(...) returns the p-value P, computed using inverse
%    interpolation into the look-up table of critical values. Small values of P
%    cast doubt on the validity of the null hypothesis. JBTEST warns when P is
%    not found within the limits of the table, i.e., outside the interval
%    [0.001, 0.50], and returns one or the other endpoint of that interval. In
%    this case, you can use the MCTOL input argument to compute a more
%    accurate value.
% 
%    [H,P,JBSTAT] = JBTEST(...) returns the test statistic JBSTAT.
% 
%    [H,P,JBSTAT,CRITVAL] = JBTEST(...) returns the critical value CRITVAL for
%    the test. When JBSTAT > CRITVAL, the null hypothesis can be rejected at a
%    significance level of ALPHA.
% 
%    [H,P,...] = JBTEST(X,ALPHA,MCTOL) computes a Monte-Carlo approximation
%    for P directly, rather than using interpolation of the pre-computed
%    tabulated values.  This is useful when ALPHA or P is outside the range of
%    the look-up table.  JBTEST chooses the number of MC replications, MCREPS,
%    large enough to make the MC standard error for P, SQRT(P*(1-P)/MCREPS),
%    less than MCTOL.
% 
%    See also LILLIETEST, KSTEST, KSTEST2, CDFPLOT.
%
%    Reference page in Doc Center
%       doc jbtest
%
%    Other functions named jbtest
%
%       ts/jbtest
%