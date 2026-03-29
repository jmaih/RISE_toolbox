% VAR Variance.
%    For vectors, Y = VAR(X) returns the variance of the values in X.  For
%    matrices, Y is a row vector containing the variance of each column of
%    X.  For N-D arrays, VAR operates along the first non-singleton
%    dimension of X.
% 
%    VAR normalizes Y by N-1 if N>1, where N is the sample size.  This is
%    an unbiased estimator of the variance of the population from which X is
%    drawn, as long as X consists of independent, identically distributed
%    samples. For N=1, Y is normalized by N. 
% 
%    Y = VAR(X,1) normalizes by N and produces the second moment of the
%    sample about its mean.  VAR(X,0) is the same as VAR(X).
% 
%    Y = VAR(X,W) computes the variance using the weight vector W.  W 
%    typically contains either counts or inverse variances.  The length of W 
%    must equal the length of the dimension over which VAR operates, and its
%    elements must be nonnegative.  If X(I) is assumed to have variance 
%    proportional to 1/W(I), then Y * MEAN(W)/W(I) is an estimate of the 
%    variance of X(I).  In other words, Y * MEAN(W) is an estimate of 
%    variance for an observation given weight 1.
% 
%    Y = VAR(X,0,"all") or Y = VAR(X,1,"all") returns the variance of all
%    elements of X. A weight of 0 normalizes by N-1 and a weight of 1 
%    normalizes by N.
% 
%    Y = VAR(X,W,DIM) takes the variance along the dimension DIM of X.
% 
%    Y = VAR(X,0,VECDIM) or Y = VAR(X,1,VECDIM) operates on the dimensions 
%    specified in the vector VECDIM. A weight of 0 normalizes by N-1 and a 
%    weight of 1 normalizes by N. For example, VAR(X,0,[1 2]) operates on
%    the elements contained in the first and second dimensions of X.
%    
%    [Y,M] = VAR(X,...) also returns the mean M of the values in X that 
%    was used to calculate the variance. If Y is the weighted variance, then
%    M is the weighted mean. 
% 
%    The variance is the square of the standard deviation (STD).
% 
%    VAR(...,NANFLAG) specifies how NaN values are treated:
% 
%    "includemissing" / "includenan"  -
%                   (default) The variance of a vector containing NaN values
%                   is also NaN.
%    "omitmissing" / "omitnan"        -
%                   Elements of X or W containing NaN values are ignored.
%                   If all elements are NaN, the result is NaN.
% 
%    Example:
%        X = [4 -2 1; 9 5 7]
%        var(X,0,1)
%        var(X,0,2)
% 
%    Class support for inputs X, W:
%       float: double, single
% 
%    See also MEAN, STD, COV, CORRCOEF.
%
%    Documentation for var
%       doc var
%
%    Other uses of var
%
%       codistributed/var
%       gpuArray/var
%       prob.BetaDistribution/var
%       prob.BinomialDistribution/var
%       prob.BirnbaumSaundersDistribution/var
%       prob.BurrDistribution/var
%       prob.ExponentialDistribution/var
%       prob.ExtremeValueDistribution/var
%       prob.GammaDistribution/var
%       prob.GeneralizedExtremeValueDistribution/var
%       prob.GeneralizedParetoDistribution/var
%       prob.HalfNormalDistribution/var
%       prob.InverseGaussianDistribution/var
%       prob.KernelDistribution/var
%       prob.LogisticDistribution/var
%       prob.LoglogisticDistribution/var
%       prob.LognormalDistribution/var
%       prob.LoguniformDistribution/var
%       prob.MultinomialDistribution/var
%       prob.NakagamiDistribution/var
%       prob.NegativeBinomialDistribution/var
%       prob.NormalDistribution/var
%       prob.PiecewiseLinearDistribution/var
%       prob.PoissonDistribution/var
%       prob.ProbabilityDistribution/var
%       prob.RayleighDistribution/var
%       prob.RicianDistribution/var
%       prob.StableDistribution/var
%       prob.tLocationScaleDistribution/var
%       prob.TriangularDistribution/var
%       prob.UniformDistribution/var
%       prob.WeibullDistribution/var
%       tabular/var
%       tall/var
%       timeseries/var
%       ts/var
%