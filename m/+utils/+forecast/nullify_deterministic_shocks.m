function shocks=nullify_deterministic_shocks(shocks,det_vars)
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

if islogical(det_vars)
    det_vars=find(det_vars);
end
if ~isempty(det_vars)
    det_vars=det_vars(:);
    nx=size(shocks,1);
    utils.forecast.check_shock_id(det_vars,nx)
    shocks(det_vars,:)=0;
end
end
