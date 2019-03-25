% &  Logical AND.
%    A & B performs a logical AND of arrays A and B and returns an array
%    containing elements set to either logical 1 (TRUE) or logical 0
%    (FALSE). An element of the output array is set to 1 if both input
%    arrays contain a non-zero element at that same array location.
%    Otherwise, that element is set to 0. A and B must have compatible
%    sizes. In the simplest cases, they can be the same size or one can be a
%    scalar. Two inputs have compatible sizes if, for every dimension, the
%    dimension sizes of the inputs are either the same or one of them is 1.
% 
%    C = AND(A,B) is called for the syntax 'A & B' when A or B is an object.
% 
%    Note that there are two logical AND operators in MATLAB. The & operator
%    performs an element-by-element AND between matrices, while the &&
%    operator performs a short-circuit AND between scalar values. See the
%    documentation for details.
% 
%    See also RELOP, OR, XOR, NOT.
%
%    Reference page in Doc Center
%       doc and
%
%    Other functions named and
%
%       codistributed/and    mtree/and     sym/and    ts/and
%       gpuArray/and         struct/and
%