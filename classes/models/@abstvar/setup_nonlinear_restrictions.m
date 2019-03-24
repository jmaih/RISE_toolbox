%--- help for generic/setup_nonlinear_restrictions ---
%
%  INTERNAL FUNCTION: Sets nonlinear restrictions
% 
%  Note:
% 
%     - uses estim_nonlinear_restrictions, which should be a cell array. Each
%     item of the array is a string of the form
%       - 'f(p1,p2,...,pn)>=h(p1,p2,...,pn)'
%       - 'f(p1,p2,...,pn)>h(p1,p2,...,pn)'
%       - 'f(p1,p2,...,pn)<=h(p1,p2,...,pn)'
%       - 'f(p1,p2,...,pn)<h(p1,p2,...,pn)'
%       - 'pj=h(p1,p2,...,pn)'
% 
%     - In the statements above,
%       - eqtn [digits|variable name]
%       - vbl [digits|variable name]
%       - lag [digits]
%       - chain [char]
%       - state [digits]
% 
%
%    Other functions named setup_nonlinear_restrictions
%
%       abstvar/setup_nonlinear_restrictions
%