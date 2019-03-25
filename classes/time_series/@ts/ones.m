% ONES   Ones array.
%    ONES(N) is an N-by-N matrix of ones.
% 
%    ONES(M,N) or ONES([M,N]) is an M-by-N matrix of ones.
% 
%    ONES(M,N,P,...) or ONES([M N P ...]) is an M-by-N-by-P-by-... array of
%    ones.
% 
%    ONES(SIZE(A)) is the same size as A and all ones.
% 
%    ONES with no arguments is the scalar 1.
% 
%    ONES(..., CLASSNAME) is an array of ones of class specified by the 
%    string CLASSNAME.
% 
%    ONES(..., 'like', Y) is an array of ones with the same data type, sparsity,
%    and complexity (real or complex) as the numeric variable Y.
% 
%    Note: The size inputs M, N, and P... should be nonnegative integers. 
%    Negative integers are treated as 0.
% 
%    Example:
%       x = ones(2,3,'int8');
% 
%    See also EYE, ZEROS.
%
%    Reference page in Doc Center
%       doc ones
%
%    Other functions named ones
%
%       codistributed/ones      codistributor2dbc/ones    gpuArray/ones
%       codistributor1d/ones    distributed/ones          ts/ones
%