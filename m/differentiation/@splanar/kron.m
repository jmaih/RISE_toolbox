% KRON   Kronecker tensor product.
%    KRON(X,Y) is the Kronecker tensor product of X and Y.
%    The result is a large matrix formed by taking all possible
%    products between the elements of X and those of Y. For
%    example, if X is 2 by 3, then KRON(X,Y) is
% 
%       [ X(1,1)*Y  X(1,2)*Y  X(1,3)*Y
%         X(2,1)*Y  X(2,2)*Y  X(2,3)*Y ]
% 
%    If either X or Y is sparse, only nonzero elements are multiplied
%    in the computation, and the result is sparse.
% 
%    Class support for inputs X,Y:
%       float: double, single
%       integers: uint8, int8, uint16, int16, uint32, int32, uint64, int64
%
%    Reference page in Doc Center
%       doc kron
%
%    Other functions named kron
%
%       splanar/kron    sym/kron
%