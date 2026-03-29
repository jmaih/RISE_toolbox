% MEAN   Average or mean value.
%    S = MEAN(X) is the mean value of the elements in X if X is a vector. 
%    For matrices, S is a row vector containing the mean value of each 
%    column. 
%    For N-D arrays, S is the mean value of the elements along the first 
%    array dimension whose size does not equal 1.
% 
%    MEAN(X,"all") is the mean of all elements in X.
% 
%    MEAN(X,DIM) takes the mean along the dimension DIM of X.
% 
%    MEAN(X,VECDIM) operates on the dimensions specified in the vector 
%    VECDIM. For example, MEAN(X,[1 2]) operates on the elements contained
%    in the first and second dimensions of X.
% 
%    S = MEAN(...,OUTTYPE) specifies the type in which the mean is performed, 
%    and the type of S. Available options are:
% 
%    "double"    -  S has class double for any input X
%    "native"    -  S has the same class as X
%    "default"   -  If X is floating point, that is double or single,
%                   S has the same class as X. If X is not floating point, 
%                   S has class double.
% 
%    S = MEAN(...,NANFLAG) specifies how NaN values are treated:
% 
%    "includemissing" / "includenan" -
%                   (default) The mean of a vector containing NaN values is NaN.
%    "omitmissing" / "omitnan"       -
%                   The mean of a vector containing NaN values is the mean
%                   of all its non-NaN elements. If all elements are NaN,
%                   the result is NaN.
% 
%    Example:
%        X = [1 2 3; 3 3 6; 4 6 8; 4 7 7]
%        mean(X,1)
%        mean(X,2)
% 
%    Class support for input X:
%       float: double, single
%       integer: uint8, int8, uint16, int16, uint32,
%                int32, uint64, int64
% 
%    See also MEDIAN, STD, MIN, MAX, VAR, COV, MODE.
%
%    Documentation for mean
%       doc mean
%
%    Other uses of mean
%
%       codistributed/mean
%       datetime/mean
%       duration/mean
%       gpuArray/mean
%       prob.BetaDistribution/mean
%       prob.BinomialDistribution/mean
%       prob.BirnbaumSaundersDistribution/mean
%       prob.BurrDistribution/mean
%       prob.ExponentialDistribution/mean
%       prob.ExtremeValueDistribution/mean
%       prob.GammaDistribution/mean
%       prob.GeneralizedExtremeValueDistribution/mean
%       prob.GeneralizedParetoDistribution/mean
%       prob.HalfNormalDistribution/mean
%       prob.InverseGaussianDistribution/mean
%       prob.KernelDistribution/mean
%       prob.LogisticDistribution/mean
%       prob.LoglogisticDistribution/mean
%       prob.LognormalDistribution/mean
%       prob.LoguniformDistribution/mean
%       prob.MultinomialDistribution/mean
%       prob.NakagamiDistribution/mean
%       prob.NegativeBinomialDistribution/mean
%       prob.NormalDistribution/mean
%       prob.PiecewiseLinearDistribution/mean
%       prob.PoissonDistribution/mean
%       prob.ProbabilityDistribution/mean
%       prob.RayleighDistribution/mean
%       prob.RicianDistribution/mean
%       prob.StableDistribution/mean
%       prob.tLocationScaleDistribution/mean
%       prob.TriangularDistribution/mean
%       prob.UniformDistribution/mean
%       prob.WeibullDistribution/mean
%       symfun/mean
%       tabular/mean
%       tall/mean
%       timeseries/mean
%       ts/mean
%