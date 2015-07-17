function y=nth_order_logistic(x,g,varargin)
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
    error('at least 3 arguments should be provided')
end
n=length(varargin);

c1=varargin{1};
x_minus_c=x-c1;
for ii=2:n
    c2=varargin{ii};
    if isa(c1,'double') && isa(c2,'double') &&  any(c2<c1)
        error('the c coefficients should be in an ascending order')
    end
    x_minus_c=x_minus_c.*(x-c2);
    c1=c2;
end
if isa(g,'double') && any(g<0)
    error('g cannot be negative')
end

y=1./(1+exp(-g.*x_minus_c));

end