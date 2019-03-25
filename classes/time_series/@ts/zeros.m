% ZEROS  Zeros array.
%    ZEROS(N) is an N-by-N matrix of zeros.
% 
%    ZEROS(M,N) or ZEROS([M,N]) is an M-by-N matrix of zeros.
% 
%    ZEROS(M,N,P,...) or ZEROS([M N P ...]) is an M-by-N-by-P-by-... array of
%    zeros.
% 
%    ZEROS(SIZE(A)) is the same size as A and all zeros.
% 
%    ZEROS with no arguments is the scalar 0.
% 
%    ZEROS(..., CLASSNAME) is an array of zeros of class specified by the
%    string CLASSNAME.
% 
%    ZEROS(..., 'like', Y) is an array of zeros with the same data type, sparsity,
%    and complexity (real or complex) as the numeric variable Y.
% 
%    Note: The size inputs M, N, and P... should be nonnegative integers. 
%    Negative integers are treated as 0.
% 
%    Example:
%       x = zeros(2,3,'int8');
% 
%    See also EYE, ONES.
%
%    Reference page in Doc Center
%       doc zeros
%
%    Other functions named zeros
%
%       codistributed/zeros      codistributor2dbc/zeros    gpuArray/zeros
%       codistributor1d/zeros    distributed/zeros          ts/zeros
%