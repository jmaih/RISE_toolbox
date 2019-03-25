%  kp_shift_solve -- Solves the problem (A1@A2@...@Ap-lambda*I)*x=b
%  where @ stands for the kronecker product
% 
%  ::
% 
% 	x = kp_shift_solve(A,lambda,b)
% 	x = kp_shift_solve(A,lambda,b,do_real)
% 
%  Args:
% 
%     - **A** [cell array]: each element Ai, i=1,2,...,p is a
%       square matrix with number of rows ni
% 
%     - **b** [vector]: right-hand side of the equality
% 
%     - **lambda** [scalar]: See formula above
% 
%     - **do_real** [true|{false}]: See formula above
% 
%  Returns:
%     :
% 
%     - **x** [N x 1 vector]: where N=n1 x n2 x ... x np
% 
%  See also : utils.kronecker.qtkp_shift_solve
%