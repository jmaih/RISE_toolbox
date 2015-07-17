function y=exponential(x,g,c)
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

if nargin<3
    c=0;
    if nargin<2
        error('at least 2 arguments should be provided')
    end
end
if isa(g,'double') && any(g<0)
    error('g cannot be negative')
end

y=1-exp(-g.*(x-c).^2);

end