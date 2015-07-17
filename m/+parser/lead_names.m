function out=lead_names(names,n)
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

if nargin<2
    n=1;
end
n=abs(n);

out=parser.concatenate_names_number(names,n);