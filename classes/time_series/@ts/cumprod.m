% CUMPROD Cumulative product of elements.
%    Y = CUMPROD(X) computes the cumulative product along the first non-singleton
%    dimension of X. Y is the same size as X.
% 
%    Y = CUMPROD(X,DIM) cumulates along the dimension specified by DIM.
%  
%    Y = CUMPROD(___,DIRECTION) cumulates in the direction specified by
%    DIRECTION using any of the above syntaxes:
% 
%      "forward" - (default) uses the forward direction, from beginning to end.
%      "reverse" -           uses the reverse direction, from end to beginning.
% 
%    Y = CUMPROD(___,NANFLAG) specifies how NaN values are treated:
% 
%    "includemissing" / "includenan" -
%                    (default) The product of elements containing NaN values
%                    is NaN.
%    "omitmissing" / "omitnan"       -
%                    The product of elements containing NaN values is the
%                    product of all non-NaN elements. If all elements are
%                    NaN, the result is 1.
% 
%    Example: 
%        X = [0 1 2; 3 4 5]
%        cumprod(X,1)
%        cumprod(X,2)
%        cumprod(X,1,"reverse")
%        cumprod(X,2,"reverse")
% 
%    See also PROD, CUMSUM, CUMMIN, CUMMAX.
%
%    Documentation for cumprod
%       doc cumprod
%
%    Other uses of cumprod
%
%       codistributed/cumprod    symbolic/cumprod    tall/cumprod
%       gpuArray/cumprod         tabular/cumprod     ts/cumprod
%       sym/cumprod
%