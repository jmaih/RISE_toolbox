function db=rand(start_date,varargin)
% Initializes a ts object with the given start data and data initialized to
%   uniform random numbers 
%
% ::
%
%   db=ts.rand(start_date,varargin)
%
% Args:
%
%    start_date (numeric | char): a valid time series (ts) date
%
%    varargin (numeric): arguments to matlab's **rand** function.
%
% Returns:
%    :
%
%    - **db** : [ts]: a time series
%
% Note:
%
%    - this is a static method and so it has to be called with the **ts.**
%      prefix
%    - ts.rand does not allow more than 3 dimensions
%
% Example:
%    ::
%
%       db=ts.rand(1990,10,1)
%       db=ts.rand('1990',10,3)
%       db=ts.rand('1990Q3',10,5,100)
%


data=rand(varargin{:});

db=ts(start_date,data);

