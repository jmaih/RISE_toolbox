function ed=distance(a,b)
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

if isempty(a)||isempty(b)
    ed=inf;
else
    dev=bsxfun(@minus,a,b);
    ed=sqrt(sum(dev.^2,1));
end
