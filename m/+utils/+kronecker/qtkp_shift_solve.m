%  qtkp_shift_solve -- Solves the problem (alpha@T1@T2@...@Tp-lambda*I)*y=c
%  where @ stands for the kronecker product for quasi-triangular matrices Ti 
% 
%  ::
% 
% 	y = qtkp_shift_solve(T,n,c,lambda)
% 	y = qtkp_shift_solve(T,n,c,lambda,alpha)
% 	y = qtkp_shift_solve(T,[],...)
% 
%  Args:
% 
%     - **T** [cell array]: each element Ti, i=1,2,...,p is a
%       quasi-triangular (square) matrix
% 
%     - **n** [vector|{[]}]: number of rows (or columns) for Ti,
%       i=1,2,...,p. 
% 
%     - **c** [vector]: right-hand side of the equality
% 
%     - **lambda** [scalar]: See formula above
% 
%     - **alpha** [scalar|2 x 2 matrix|{1}]: First element of the kronecker
%       product
% 
%  Returns:
%     :
% 
%     - **y** [N x 1 vector]: where N=n1 x n2 x ... x np
% 
%  See also : utils.kronecker.kp_shift_solve
%