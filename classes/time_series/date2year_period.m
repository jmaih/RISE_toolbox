function [year,period,freq,frequency]=date2year_period(date)
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

[dec,success]=decompose_date(date);

if ~success
    
    error('wrong date')
    
end

year=[dec.year];

year=year(:);

period=[dec.period];

period=period(:);

freq=dec(1).freq;

frequency=dec(1).frequency;

end