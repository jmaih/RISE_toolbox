function Y=range(this,varargin)
% Vverloads Matlab's RANGE for ts objects. returns the range of the
% values in the time series
%
% ::
%
%    Y=RANGE(this,varargin)
%
% Args:
%
%    - **this** [ts]: time series object
%    - **varargin**: additional matlab arguments for the RANGE function
%
% Returns:
%
%    - **Y** [numeric]: Difference between maximum and minimum values
%
% See also:
%    - range

Y=range(this.data,varargin{:});

end