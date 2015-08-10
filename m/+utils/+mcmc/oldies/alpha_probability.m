function a=alpha_probability(f_theta,f0)
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

r=exp(f_theta-f0);
a=min(1,r);