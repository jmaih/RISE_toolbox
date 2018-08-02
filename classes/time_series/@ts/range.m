function Y=range(this,varargin)
% Overloads Matlab's range for ts objects. returns the range of the
% values in the time series
%
% ::
%
%    Y = range(this,varargin);
%
% Args:
%
%    this (ts object): time series object
%    varargin : additional matlab arguments for the range function
%
% Returns:
%    :
%
%    - **Y** [numeric]: Difference between maximum and minimum values
%
% See also:
%    - range
%

Y=range(this.data,varargin{:});

end