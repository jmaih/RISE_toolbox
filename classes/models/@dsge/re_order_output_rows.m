function y=re_order_output_rows(obj,y)
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


y=reshape(y(obj.inv_order_var,:,:),size(y));

end