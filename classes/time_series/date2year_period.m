function [year,period,freq,frequency]=date2year_period(date)
% H1 line
%
% Syntax
% -------
% ::
%
% Inputs
% -------
%
% Outputs
% --------
%
% More About
% ------------
%
% Examples
% ---------
%
% See also:

[dec,success]=decompose_date(date);

if ~success
    
    error('wrong date')
    
end

year=cellfun(@str2double,{dec.year},'uniformOutput',true);

year=year(:);

period=cellfun(@str2double,{dec.period},'uniformOutput',true);

period=period(:);

freq=dec(1).freq;

frequency=dec(1).frequency;

end