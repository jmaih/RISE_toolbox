function db=ones(start_date,varargin)
% ones overloads ones for ts objects
%
% ::
%
%   db=ts.ones(start_date,varargin)
%
% Args:
%
%    - **start_date** : [numeric|char]: a valid time series (ts) date
%
%    - varargin : [numeric]: arguments to matlab's **ones** function.
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
%
%    - ts.ones does not allow more than 3 dimensions
%
% Example:
%
%      db=ts.ones(1990,10,1)
%      db=ts.ones('1990',10,3)
%      db=ts.ones('1990Q3',10,5,100)
%
%    See also:


data=ones(varargin{:});

db=ts(start_date,data);

