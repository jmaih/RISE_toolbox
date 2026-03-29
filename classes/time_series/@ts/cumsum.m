% CUMSUM Cumulative sum of elements.
%    Y = CUMSUM(X) computes the cumulative sum along the first non-singleton
%    dimension of X. Y is the same size as X.
%  
%    Y = CUMSUM(X,DIM) cumulates along the dimension specified by DIM.
%  
%    Y = CUMSUM(___,DIRECTION) cumulates in the direction specified by
%    DIRECTION using any of the above syntaxes:
% 
%      "forward" - (default) uses the forward direction, from beginning to end.
%      "reverse" -           uses the reverse direction, from end to beginning.
% 
%    Y = CUMSUM(___,NANFLAG) specifies how NaN values are treated:
% 
%    "includemissing" / "includenan" -
%                    (default) The sum of elements containing NaN values is NaN.
%    "omitmissing" / "omitnan"       -
%                    The sum of elements containing NaN values is the sum of
%                    all non-NaN elements. If all elements are NaN, the result
%                    is 0.
% 
%    Example: 
%        X = [0 1 2; 3 4 5]
%        cumsum(X,1)
%        cumsum(X,2)
%        cumsum(X,1,"reverse")
%        cumsum(X,2,"reverse")
% 
%    See also SUM, MOVSUM, CUMPROD, CUMMIN, CUMMAX.
%
%    Documentation for cumsum
%       doc cumsum
%
%    Other uses of cumsum
%
%       codistributed/cumsum    sym/cumsum         tall/cumsum
%       duration/cumsum         symbolic/cumsum    ts/cumsum
%       gpuArray/cumsum         tabular/cumsum
%