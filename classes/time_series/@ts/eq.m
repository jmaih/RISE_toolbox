% ==  Equal.
%    A == B does element by element comparisons between A and B and returns
%    an array with elements set to logical 1 (TRUE) where the relation is
%    true and elements set to logical 0 (FALSE) where it is not. A and B
%    must have compatible sizes. In the simplest cases, they can be the same
%    size or one can be a scalar. Two inputs have compatible sizes if, for
%    every dimension, the dimension sizes of the inputs are either the same
%    or one of them is 1.
% 
%    C = EQ(A,B) is called for the syntax 'A == B' when A or B is an object.
%
%    Reference page in Doc Center
%       doc eq
%
%    Other functions named eq
%
%       categorical/eq      MException/eq      string/eq
%       codistributed/eq    mtree/eq           sym/eq
%       fints/eq            opaque/eq          timer/eq
%       gpuArray/eq         scribehandle/eq    timeseries/eq
%       handle/eq           serial/eq          ts/eq
%       instrument/eq
%