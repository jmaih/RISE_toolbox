function y=re_order_output_rows(obj,y)
% INTERNAL FUNCTION
%

y=reshape(y(obj.inv_order_var,:,:),size(y));

end