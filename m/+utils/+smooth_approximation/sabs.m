function a=sabs(x,tol)
% H1 line
%
% ::
%
%
% Args:
%
% Returns:
%    :
%
% Note:
%
% Example:
%
%    See also:

if nargin<2
    tol=1e-5;
end
% smooth approximation of abs
a=x.^2./sqrt(x.^2+tol.^2);
end