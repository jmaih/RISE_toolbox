% CORR Linear or rank correlation.
%    RHO = CORR(X) returns a P-by-P matrix containing the pairwise linear
%    correlation coefficient between each pair of columns in the N-by-P
%    matrix X.
% 
%    RHO = CORR(X,Y,...) returns a P1-by-P2 matrix containing the pairwise
%    correlation coefficient between each pair of columns in the N-by-P1 and
%    N-by-P2 matrices X and Y.
% 
%    [RHO,PVAL] = CORR(...) also returns PVAL, a matrix of p-values for
%    testing the hypothesis of no correlation against the alternative that
%    there is a non-zero correlation.  Each element of PVAL is the p-value
%    for the corresponding element of RHO.  If PVAL(i,j) is small, say less
%    than 0.05, then the correlation RHO(i,j) is significantly different
%    from zero.
% 
%    [...] = CORR(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies additional
%    parameters and their values.  Valid parameters are the following:
% 
%         Parameter  Value
%          'type'    'Pearson' (the default) to compute Pearson's linear
%                    correlation coefficient, 'Kendall' to compute Kendall's
%                    tau, or 'Spearman' to compute Spearman's rho.
%          'rows'    'all' (default) to use all rows regardless of missing
%                    values (NaNs), 'complete' to use only rows with no
%                    missing values, or 'pairwise' to compute RHO(i,j) using
%                    rows with no missing values in column i or j.
%          'tail'    The alternative hypothesis against which to compute
%                    p-values for testing the hypothesis of no correlation.
%                    Choices are:
%                       TAIL         Alternative Hypothesis
%                    ---------------------------------------------------
%                      'both'     correlation is not zero (the default)
%                      'right'    correlation is greater than zero
%                      'left'     correlation is less than zero
% 
%    The 'pairwise' option for the 'rows' parameter can produce RHO that is
%    not positive semi-definite.  The 'complete' option always produces a
%    positive semi-definite RHO, but when data are missing, the estimates
%    may be based on fewer observations.
% 
%    CORR computes p-values for Pearson's correlation using a Student's t
%    distribution for a transformation of the correlation.  This is exact
%    when X and Y are normal.  CORR computes p-values for Kendall's tau and
%    Spearman's rho using either the exact permutation distributions (for
%    small sample sizes), or large-sample approximations.
% 
%    CORR computes p-values for the two-tailed test by doubling the more
%    significant of the two one-tailed p-values.
% 
%    See also CORRCOEF, PARTIALCORR, TIEDRANK.
%
%    Documentation for corr
%       doc corr
%
%    Other uses of corr
%
%       gpuArray/corr    ssm/corr    tall/corr    ts/corr
%