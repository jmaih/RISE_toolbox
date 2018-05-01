function a=alpha_probability(f_theta,f0)
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

r=exp(f_theta-f0);
a=min(1,r);