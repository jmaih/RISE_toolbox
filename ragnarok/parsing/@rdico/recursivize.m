%--- help for rdico.recursivize ---
%
% RECURSIVIZE Ensures all operations in the cell array are recursive.
%    operations = RECURSIVIZE(operations) takes a cell array of operations 
%    and ensures that each operation is recursive. If an operation is 
%    non-recursive, it replaces the LHS of the non-recursive operation with
%    its RHS in all the previous operations where it appears.
% 
%    Input:
%    operations - Cell array of strings, where each string is an equation
%                 in the format 'lhs = rhs'.
% 
%    Output:
%    operations - Cell array of strings, with non-recursive operations
%                 modified to be recursive by substituting their LHS with
%                 their RHS in the previous equations.
%  
%  Example usage
%  operations = {
%      'x = y + z';    % non-recursive
%      'w = u + v';    % recursive
%      'y = u + w';    % recursive
%      'z = x + 1';    % non-recursive
%  };
%  
%  final_operations = recursivize(operations);
%