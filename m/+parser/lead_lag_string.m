function add_str=lead_lag_string(lead_or_lag)
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

if lead_or_lag>0
    add_str='LEAD_';
elseif lead_or_lag<0
    add_str='LAG_';
else
    error('input must be strictly positive or strictly negative')
end
end