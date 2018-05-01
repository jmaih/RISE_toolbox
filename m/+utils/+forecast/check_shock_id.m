function check_shock_id(id,nx)
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

id=id(:);
if any(id>nx)||...
        any(id<=0)||...
        any(floor(id)~=id)
    msg=['id should be an integer between 1 and ',int2str(nx)];
    error(msg)
end
end
