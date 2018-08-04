function Y=range(this,varargin)
% RANGE overloads Matlab's RANGE for ts objects. returns the range of the
% values in the time series
%
% ::
%
%    Y=RANGE(this,varargin)
%
% Args:
%
%    - **this** [ts]: time series object
%
%    - **varargin**: additional matlab arguments for the RANGE function
%
% Output Args:
%
%    - **Y** [numeric]: Difference between maximum and minimum values
%
% See also : RANGE

Y=range(this.data,varargin{:});

end