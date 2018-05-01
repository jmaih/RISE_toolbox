function Y=range(this,varargin)
% H1 line
%
% ::
%
%
% Args:
%
% Returns:
%    :
%
% Note:
%
% Example:
%
%    See also:

% range  Sample range.
%    Y = range(db) returns the range of the values in db.  For a vector input,
%    Y is the difference between the maximum and minimum values.  For a
%    matrix input, Y is a vector containing the range for each column.  For
%    N-D arrays, range operates along the first non-singleton dimension.
% 
%    range treats NaNs as missing values, and ignores them.
% 
%    Y = range(db,DIM) operates along the dimension DIM.
% 
%    See also iqr, mad, max, min, std.

Y=range(this.data,varargin{:});

end