function y = nansum(x,dim)
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

x(isnan(x)) = 0;
if nargin == 1 
    y = sum(x);
else 
    y = sum(x,dim);
end
