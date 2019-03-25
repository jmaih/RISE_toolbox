%  INTERNAL FUNCTION: creates blocks of equation-variables so as to permit a recursive solution
% 
%  ::
% 
%    [equation_blocks,variable_blocks]=block_triangularize(Incidence,do_blocks)
% 
%  Args:
% 
%     - **Incidence** [sparse]: Incidence matrix where rows represent equations
%       and columns represent variables.
%     - **do_blocks** [false|{true}]: triggers the computation of blocks.
% 
%  Returns:
%     :
% 
%     - **equation_blocks** [cell array]: blocks of equations IDs
%     - **variable_blocks** [cell array]: blocks of variables IDs
% 
%