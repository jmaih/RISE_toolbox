%  INTERNAL FUNCTION: Computes the product between a matrix and a sum of tensor products
% 
%  ::
% 
%    res=A_times_sum_perms(A,kron_matrices,matsizes,add_lead_term,P1,...,Pm)
% 
%  Args:
% 
%     - **A** [matrix]:
%     - **kron_matrices** [cell array]: e.g. {M1,M2,...,Mn}, where the Mi's are
%       the matrices entering the various tensor products.
%     - **matsizes** [k x 2 vector]: sizes of the matrices. Note need to be
%       taken of the fact that these sizes are not necessarily the sizes of the
%       matrices entering **kron_matrices** . e.g. suppose we want to compute
%       A(UUB+UBU+BUU), where U is unknown but UU is known. Then we can have
%       kron_matrices={UU,B} but matsizes=[size(U);size(U);size(B)]
%     - **add_lead_term** [{true}|false]: if true, the term A(M1*M2*...*Mn) is
%       added to the sum
%     - **P1,...,Pm** [vectors]: permutations of the tensor products entering
%       the sum. The number of elements entering each Pi has to be the same as
%       the number of rows of **matsizes**, which represents the effective number
%       of matrices.
% 
%  Returns:
%     :
% 
%     - **res** [matrix]: result of A(B*C*D+C*D*B+...)
% 
%