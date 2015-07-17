function a=smin(x,tol)
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

% smooth approximation of min(x,0)
if nargin<2
    tol=1e-5;
end
a=0.5*(x-utils.smooth_approximation.sabs(x,tol));
end