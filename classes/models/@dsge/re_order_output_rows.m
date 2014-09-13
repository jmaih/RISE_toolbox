function y=re_order_output_rows(obj,y)

y=reshape(y(obj.inv_order_var.after_solve,:),size(y));

end