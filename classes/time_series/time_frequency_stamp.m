function [stamp,unstamp]=time_frequency_stamp()
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

% stamps helps put a stamp on the serial numbers such that the frequency is
% recoverable

% unstamp inverts the stamp above

option=2;
switch option
    case 1
        stamp=@(x)0.01*x;
        unstamp=@(x)round(100*x);
    case 2
        stamp=@(x)1./(1+x);
        unstamp=@(x)round(1./x-1);
    otherwise
        error('case not implemented')
end


