function vcov=farthest(vcov0,debug)

% farthest -- computes farthest covariance matrix when failure
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
    debug=false;
end

vcov=utils.cov.nearest(vcov0,debug,true);

end