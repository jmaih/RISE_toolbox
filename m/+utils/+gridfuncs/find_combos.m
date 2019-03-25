%  INTERNAL FUNCTION: find unique permutations arising from tensors of sums
% 
%  ::
% 
%    c=find_combos(acell,bcell)
% 
%  Args:
% 
%     - **acell** [2x1 or 1x2 cell]: where the first element is a number or a
%       char and the second element is an integer describing the number of times
%       the element in the first cell occurs
%     - **bcell** [2x1 or 1x2 cell]: where the first element is a number or a
%       char and the second element is an integer describing the number of times
%       the element in the first cell occurs
% 
%  Returns:
%     :
% 
%     - **c** [n x k matrix]: char or doubles depending on the type of inputs
% 
%  Note:
% 
%     - using the Pascal's triangle, we know that
%       (a+b)^5=a^5+5a^4*b+10a^3*b^2+10a^2*b^3+5a*b^4+b^5. Suppose now that a and
%       b are matrices and let's focus on one term, say 10a^3*b^2. We know that
%       we will have 10 combinations of kronecker products in which a appears 3
%       times and b appears twice. But then how to find them? This routine is
%       designed to solve that problem efficiently.
% 
%  Example:
% 
%     - suppose we want to find all permutations in which u appears 4 times and
%       h appears twice. c=find_combos({'u',4},{'h',2}),size(c)
%     - suppose we want to find all permutations in which 1 appears 4 times and
%       2 appears three times. c=find_combos({1,4},{2,3})
% 
%