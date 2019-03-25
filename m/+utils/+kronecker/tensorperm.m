%  INTERNAL FUNCTION: finds indices (rows and columns) for the permutations of a
%  tensor
% 
%  ::
% 
%    [irows,jcols]=tensorperm(matsizes,order_1,order_2,...,order_n)
%    [irows,jcols]=tensorperm(matsizes,order_1,order_2,...,order_n,'grid')
% 
%  Args:
% 
%     - **matsizes** [matrix]: k x 2 matrix of size of the various matrices
%       entering the tensor. Each row represents the size of a matrix and it is
%       assumed that the main(or first) tensor product is ordered [1,2,...,k]
%     - **order_i** [vector]: permutation of [1,2,...,k]. N.B: all elements
%       1,2,...,k should be part of the vector and the vector should have exactly
%       k elements
%     - **grid** [string]: if present, a grid is used to compute the indexes of
%       the main kronecker product
% 
%  Returns:
%     :
% 
%     - **irows** [vector|matrix]: each column corresponds to the rows for a
%       particular ordering as listed in varargin
%     - **jcols** [vector|matrix]: each column corresponds to the columns for a
%       particular ordering as listed in varargin
% 
%  See also:
%     sum_permutations
% 
%