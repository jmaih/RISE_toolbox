function y=re_order_output_rows(obj,y)
% H1 line
%
% Syntax
% -------
% ::
%
% Inputs
% -------
%
% Outputs
% --------
%
% More About
% ------------
%
% Examples
% ---------
%
% See also: 


y=reshape(y(obj.inv_order_var.after_solve,:),size(y));

end