%--- help for sym/collect ---
%
% COLLECT Collect coefficients.
%    COLLECT(S,v) regards each element of the symbolic matrix S as a
%    polynomial in v and rewrites S in terms of the powers of v.
%    COLLECT(S) uses the default variable determined by SYMVAR.
%    COLLECT(S, 'f') uses all symbolic calls to the function f as variables.
%    COLLECT(S, {'f1', ..., 'fk'}) uses all symbolic calls to any of the
%    functions f1, ..., fk as variables.
% 
%    Examples:
%       syms x y
% 
%       collect(x^2*y + y*x - x^2 - 2*x)  returns (y - 1)*x^2 + (y - 2)*x
% 
%       f = -1/4*x*exp(-2*x)+3/16*exp(-2*x)
%       collect(f,exp(-2*x))  returns -(x/4 - 3/16)/exp(2*x)
% 
%       f = x*sin(x) + sin(x) + x*sin(2*x) - sin(2*x)
%       collect(f, 'sin')  returns (x + 1)*sin(x) + (x - 1)*sin(2*x)
% 
%    See also SYM/SIMPLIFY, SYM/FACTOR, SYM/EXPAND, SYM/SYMVAR.
%
%    Other functions named collect
%
%       ts/collect
%