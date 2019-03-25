% REGRESS Multiple linear regression using least squares.
%    B = REGRESS(Y,X) returns the vector B of regression coefficients in the
%    linear model Y = X*B.  X is an n-by-p design matrix, with rows
%    corresponding to observations and columns to predictor variables.  Y is
%    an n-by-1 vector of response observations.
% 
%    [B,BINT] = REGRESS(Y,X) returns a matrix BINT of 95% confidence
%    intervals for B.
% 
%    [B,BINT,R] = REGRESS(Y,X) returns a vector R of residuals.
% 
%    [B,BINT,R,RINT] = REGRESS(Y,X) returns a matrix RINT of intervals that
%    can be used to diagnose outliers.  If RINT(i,:) does not contain zero,
%    then the i-th residual is larger than would be expected, at the 5%
%    significance level.  This is evidence that the I-th observation is an
%    outlier.
% 
%    [B,BINT,R,RINT,STATS] = REGRESS(Y,X) returns a vector STATS containing, in
%    the following order, the R-square statistic, the F statistic and p value
%    for the full model, and an estimate of the error variance.
% 
%    [...] = REGRESS(Y,X,ALPHA) uses a 100*(1-ALPHA)% confidence level to
%    compute BINT, and a (100*ALPHA)% significance level to compute RINT.
% 
%    X should include a column of ones so that the model contains a constant
%    term.  The F statistic and p value are computed under the assumption
%    that the model contains a constant term, and they are not correct for
%    models without a constant.  The R-square value is one minus the ratio of
%    the error sum of squares to the total sum of squares.  This value can
%    be negative for models without a constant, which indicates that the
%    model is not appropriate for the data.
% 
%    If columns of X are linearly dependent, REGRESS sets the maximum
%    possible number of elements of B to zero to obtain a "basic solution",
%    and returns zeros in elements of BINT corresponding to the zero
%    elements of B.
% 
%    REGRESS treats NaNs in X or Y as missing values, and removes them.
% 
%    See also LSCOV, POLYFIT, REGSTATS, ROBUSTFIT, STEPWISE.
%
%    Reference page in Doc Center
%       doc regress
%
%    Other functions named regress
%
%       ts/regress
%