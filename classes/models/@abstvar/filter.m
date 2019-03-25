% FILTER One-dimensional digital filter.
%    Y = FILTER(B,A,X) filters the data in vector X with the
%    filter described by vectors A and B to create the filtered
%    data Y.  The filter is a "Direct Form II Transposed"
%    implementation of the standard difference equation:
% 
%    a(1)*y(n) = b(1)*x(n) + b(2)*x(n-1) + ... + b(nb+1)*x(n-nb)
%                          - a(2)*y(n-1) - ... - a(na+1)*y(n-na)
% 
%    If a(1) is not equal to 1, FILTER normalizes the filter
%    coefficients by a(1). 
% 
%    FILTER always operates along the first non-singleton dimension,
%    namely dimension 1 for column vectors and non-trivial matrices,
%    and dimension 2 for row vectors.
% 
%    [Y,Zf] = FILTER(B,A,X,Zi) gives access to initial and final
%    conditions, Zi and Zf, of the delays.  Zi is a vector of length
%    MAX(LENGTH(A),LENGTH(B))-1, or an array with the leading dimension 
%    of size MAX(LENGTH(A),LENGTH(B))-1 and with remaining dimensions 
%    matching those of X.
% 
%    FILTER(B,A,X,[],DIM) or FILTER(B,A,X,Zi,DIM) operates along the
%    dimension DIM.
% 
%    Tip:  If you have the Signal Processing Toolbox, you can design a
%    filter, D, using DESIGNFILT.  Then you can use Y = FILTER(D,X) to
%    filter your data.
% 
%    See also FILTER2, FILTFILT, FILTIC, DESIGNFILT.
% 
%    Note: FILTFILT, FILTIC and DESIGNFILT are in the Signal Processing
%    Toolbox.
%
%    Reference page in Doc Center
%       doc filter
%
%    Other functions named filter
%
%       abstvar/filter    garch/filter       statespace/filter
%       arima/filter      gjr/filter         tall/filter
%       dsge/filter       gpuArray/filter    timeseries/filter
%       dssm/filter       LagOp/filter       varm/filter
%       egarch/filter     regARIMA/filter    vecm/filter
%       fints/filter      ssm/filter
%