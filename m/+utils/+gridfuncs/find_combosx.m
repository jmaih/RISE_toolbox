%  find_combosx -- find unique permutations arising from tensors of sums
% 
%  ::
% 
% 
%    c=find_combosx(a1_cell,a2_cell,...,an_cell)
% 
%  Args:
% 
%     - **ai_cell** [1x2 cell]: where the first element is a number or a
%     char and the second element is an integer describing the number of times
%     the element in the first cell occurs
% 
%  Returns:
%     :
% 
%     - **c** [n x k matrix]: char or doubles depending on the type of inputs
% 
%     - **p** [n x k matrix]: matrix of permutations
% 
%     - **pstar** [n x k matrix]: matrix of permutations that can be used in
%     the reordering of kronecker products. It is p above but without the
%     multiplicity
% 
%  Note:
% 
%     - using the Pascal's triangle, we know that
%     (a+b)^5=a^5+5a^4*b+10a^3*b^2+10a^2*b^3+5a*b^4+b^5. Suppose now that a and
%     b are matrices and let's focus on one term, say 10a^3*b^2. We know that
%     we will have 10 combinations of kronecker products in which a appears 3
%     times and b appears twice. But then how to find them? This routine is
%     designed to solve that problem efficiently.
% 
%  Example:
% 
%     - suppose we want to find all permutations in which Q appears 3 times,
%       P appears 2 times and O appears 2 times.
%       c=find_combosx({'Q',3},{'P',2},{'O',2}),size(c) 
% 
%     - suppose we want to find all permutations in which 5 appears 3 times,
%       7 appears 2 times and 9 appears 2 times.
%       c=find_combosx({5,3},{7,2},{9,2}),size(c) 
% 
%     See also: find_combos
%