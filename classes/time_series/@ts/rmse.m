% RMSE Root-mean-square error.
%    E = RMSE(F,A) computes the root-mean-square error between arrays F and
%    A. F is the forecast (predicted) data, and A is the actual (observed)
%    data. F and A must be double or single arrays.
% 
%    E = RMSE(F,A,"all") computes the root-mean-square error of all elements
%    in F and A.
% 
%    E = RMSE(F,A,DIM) operates along the dimension DIM of F and A.
% 
%    E = RMSE(F,A,VECDIM) operates on the dimensions specified in the vector
%    VECDIM. For example, RMSE(F,A,[1,2]) operates on the elements contained
%    in the first and second dimensions of F and A.
% 
%    E = RMSE(...,NANFLAG) specifies how NaN values are treated:
% 
%    "includemissing" / "includenan" - 
%                   (default) If F or A contains NaN values, the result is 
%                   also NaN.
%    "omitmissing" / "omitnan"       -
%                   Elements of F, A, or weights W containing NaN values are
%                   ignored. If all elements of F, A, or W are NaN, the
%                   result is NaN.
% 
%    E = RMSE(...,'Weight',W) computes the weighted error between F and A
%    using the array of weights W. The elements of W must be nonnegative. If
%    W is a vector, its length must equal the length of the dimension over
%    which RMSE is operating. If W is a matrix or multidimensional array, it
%    must have the same dimensions as F, A, or F-A.  The weight argument is
%    not supported with a vector dimension argument or the "all" flag.
% 
%    When F or A is complex, the root-mean-square error is computed using
%    the complex magnitude of F-A.
% 
%    Example:
%        F = [1 2 2;10 5 10;9 10 7]; % three forecasts
%        A = [1;9;10]; % actual
%        E = rmse(F,A)
% 
%    See also MAPE, MEAN, RMS.
%
%    Documentation for rmse
%       doc rmse
%
%    Other uses of rmse
%
%       distributed/rmse    tall/rmse    ts/rmse
%