function [flag,year,period,freq,frequency]=is_serial(x)
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


[freq,frequency,year,period]=serial2frequency(x);

flag=~isempty(freq);

end