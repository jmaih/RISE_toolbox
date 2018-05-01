function shocks=replace_impulse(shocks,shock_id,k_plus_1,new_impulse)
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


nx=size(shocks,1);
utils.forecast.check_shock_id(shock_id,nx);

shocks(shock_id,k_plus_1)=new_impulse;

end