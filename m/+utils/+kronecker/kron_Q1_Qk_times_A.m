%  INTERNAL FUNCTION: Computes (Q1*Q2*...*Qk)A where * denotes the kronecker product
% 
%  ::
% 
%     C=kron_Q1_Qk_times_A(A,Q1,Q2,...,Qk)
% 
%  Args:
% 
%     - **A** [matrix]: matrix on the left-hand side
%     - **Qi** [matrix]: matrix on the kronecker block
% 
%  Returns:
%     :
% 
%     - **C** [matrix]: result
% 
%  See also:
%     A_times_kron_Q1_Qk
% 
%