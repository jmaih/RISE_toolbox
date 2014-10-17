function ed=distance(a,b)
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

if isempty(a)||isempty(b)
    ed=inf;
else
    dev=bsxfun(@minus,a,b);
    ed=sqrt(sum(dev.^2,1));
end
