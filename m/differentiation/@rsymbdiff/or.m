% |   Logical OR.
%    A | B performs a logical OR of arrays A and B and returns an array
%    containing elements set to either logical 1 (TRUE) or logical 0
%    (FALSE). An element of the output array is set to 1 if either input
%    array contains a non-zero element at that same array location.
%    Otherwise, that element is set to 0. A and B must have compatible
%    sizes. In the simplest cases, they can be the same size or one can be a
%    scalar. Two inputs have compatible sizes if, for every dimension, the
%    dimension sizes of the inputs are either the same or one of them is 1.
% 
%    C = OR(A,B) is called for the syntax 'A | B' when A or B is an object.
% 
%    Note that there are two logical OR operators in MATLAB. The | operator
%    performs an element-by-element OR between matrices, while the ||
%    operator performs a short-circuit OR between scalar values.
% 
%    See <a href="matlab:helpview('matlab','MATLAB_OPS')">MATLAB Operators and Special Characters</a> for more details.
% 
%    See also AND, XOR, NOT.
%
%    Documentation for or
%       doc or
%
%    Other uses of or
%
%       adolm/or
%       aplanar/or
%       codistributed/or
%       gpuArray/or
%       matlab.unittest.constraints.AbsoluteTolerance/or
%       matlab.unittest.constraints.BooleanConstraint/or
%       matlab.unittest.constraints.StartsWithSubstring/or
%       matlab.unittest.internal.selectors.Modifier/or
%       matlab.unittest.internal.selectors.NeverFilterSelector/or
%       mtree/or
%       rsymbdiff/or
%       splanar/or
%       sym/or
%       symbolic/or
%       tabular/or
%