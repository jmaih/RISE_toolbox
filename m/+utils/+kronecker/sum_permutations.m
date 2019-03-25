%  INTERNAL FUNCTION: sums the permutations of a tensor product
% 
%  ::
% 
%    P=sum_permutations(ABCD,matsizes,[],order_1,order_2,...,order_n)
%    P=sum_permutations(ABCD,matsizes,[],{order_1,order_2,...,order_n})
%    P=sum_permutations(ABCD,matsizes,options,order_1,order_2,...,order_n)
%    P=sum_permutations(ABCD,matsizes,options,{order_1,order_2,...,order_n})
% 
%  Args:
% 
%     - **ABCD** [matrix]: tensor product of k matrices
%     - **matsizes** [matrix]: k x 2 matrix of size of the various matrices
%       entering the tensor. Each row represents the size of a matrix and it is
%       assumed that the main(or first) tensor product is ordered [1,2,...,k]
%     - **options** [empty|struct]: structure with various options such as
% 
%       - **algo** [{'shuf1'}|'shuf2'|'old'|'perm']: shuf1 shuffles the
%         matrix, shuf2 pre and post multiplies the matrix, old may or may not
%         construct a grid, perm uses a permutation-type of strategy.
%       - **use_grid** [true|{false}]: use grid in the old algorithm:a grid is
%         used to compute the indexes of the main kronecker product
%       - **skip_first** [true|{false}]: if true, the original input matrix is
%         not added to the sum
% 
%     - **order_i** [vector]: permutation of [1,2,...,k]. N.B: all elements
%       1,2,...,k should be part of the vector and the vector should have exactly
%       k elements
% 
%  Returns:
%     :
% 
%     - **P** sum of ABCD and its permutations
% 
%  See also:
%     tensorperm
% 
%