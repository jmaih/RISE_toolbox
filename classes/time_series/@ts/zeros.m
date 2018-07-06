function db=zeros(start_date,varargin)
% Initializes a ts object with the given start data and data initialized to zeros
%
% ::
%
%   db=ts.zeros(start_date,varargin)
%
% Args:
%
%    start_date (numeric | char): a valid time series (ts) date
%    varargin (numeric): arguments to matlab's **zeros** function.
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
%    - ts.zeros does not allow more than 3 dimensions
%
% Example:
%    ::
%
%       db=ts.zeros(1990,10,1)
%       db=ts.zeros('1990',10,3)
%       db=ts.zeros('1990Q3',10,5,100)
%

data=zeros(varargin{:});

db=ts(start_date,data);
