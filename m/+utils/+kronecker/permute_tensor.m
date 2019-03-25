%  INTERNAL FUNCTION: Permute a kronecker product
% 
%  ::
% 
%    [H1,H2,...,Hn]=permute_tensor(A1_Ak,matsizes,order_1,...,order_n)
% 
%  Args:
% 
%     - **A1_Ak** [matrix]: representing the kronecker product A1*A2*A3*...*Ak
%     - **matsizes** [matrix]: k x 2 matrix of size of the various matrices
%       entering the tensor. Each row represents the size of a matrix and it is
%       assumed that the main(or first) tensor product is ordered [1,2,...,k]
%     - **order_i** [vector]: permutation of [1,2,...,k]. N.B: all elements
%       1,2,...,k should be part of the vector and the vector should have exactly
%       k elements
% 
%  Returns:
%     :
% 
%     - **Hi** [matrix]: tensor permutation of A1_Ak according to order "i"
% 
%  See also:
%     - sum_permutations
%     - tensorperm
% 
%