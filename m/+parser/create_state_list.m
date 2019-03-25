%--- help for dsge/create_state_list ---
%
%  INTERNAL FUNCTION: create_state_list creates the list of the state variables in the solution
% 
%  ::
% 
%    final_list=create_state_list(m)
%    final_list=create_state_list(m,orders)
% 
%  Args:
% 
%     m (dsge | rise): model object
% 
%     orders (integer array | {1:m.options.solve_order}): approximation
%       orders
% 
%     compact_form (true | {false}): if true, only unique combinations
%       will be returned. Else, all combinations will be returned.
% 
%  Returns:
%     :
% 
%     - **final_list** [cellstr] : list of the state variables
%     - **kept** [vector] : location of kept state variables (computed only if
%       compact_form is set to true)
% 
%