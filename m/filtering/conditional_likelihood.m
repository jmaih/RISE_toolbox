function [lik,v_iF_v]=conditional_likelihood(v,iF,dF,pp)
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

v_iF_v=v'*iF*v;
lik=(2*pi)^(-.5*pp)*dF^(-.5)*exp(-.5*v_iF_v);
end
