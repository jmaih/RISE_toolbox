%--- help for dsge/state_var_list ---
%
%  INTERNAL FUNCTION: state_var_list creates the list of the state variables in the solution
% 
%  ::
% 
%    final_list=state_var_list(m)
%    final_list=state_var_list(m,orders)
%    final_list=state_var_list(m,orders,compact_form)
% 
%  Args:
% 
%     - m (dsge | rise): model object
% 
%     - orders (integer array | {1:m.options.solve_order}): approximation
%       orders
% 
%     - compact_form (true | {false}): if true, only unique combinations
%       will be returned. Else, all combinations will be returned.
% 
%  Returns:
%     :
% 
%     - **final_list** [cellstr] : list of the state variables
%     - **kept** [vector] : location of kept state variables (computed only if
%       compact_form is set to true)
%