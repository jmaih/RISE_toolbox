function m = nanmean(x,dim)
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

nans = isnan(x);
x(nans) = 0;

if nargin == 1 
    n = sum(~nans);
    n(n==0) = NaN;
    m = sum(x) ./ n;
else
    n = sum(~nans,dim);
    n(n==0) = NaN; 
    m = sum(x,dim) ./ n;
end

