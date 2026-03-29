% RESHAPE Reshape array by rearranging existing elements
%    B = RESHAPE(A,M,N) or RESHAPE(A,[M,N]) returns the M-by-N matrix whose
%    elements are taken columnwise from A. A must have M*N elements.
% 
%    B = RESHAPE(A,M,N,P,...) or RESHAPE(A,[M,N,P,...]) returns a
%    multidimensional array with the same elements as A but reshaped to have
%    the size M-by-N-by-P-by-.... The product of the specified dimensions,
%    M*N*P*..., must be the same as NUMEL(A).
% 
%    B = RESHAPE(A,...,[],...) reshapes A according to the specified
%    dimension sizes. You can specify a single dimension size of [] to have
%    the dimension size automatically calculated, such that the product of
%    the dimensions equals NUMEL(A). The value of NUMEL(A) must be evenly
%    divisible by the product of the specified dimensions.
% 
%    See also SQUEEZE, SHIFTDIM, COLON, RESIZE.
%
%    Documentation for reshape
%       doc reshape
%
%    Other uses of reshape
%
%       calendarDuration/reshape
%       categorical/reshape
%       codistributed/reshape
%       datetime/reshape
%       digraph/reshape
%       duration/reshape
%       gpuArray/reshape
%       lib.pointer/reshape
%       matlab.mixin.indexing.RedefinesParen/reshape
%       RandStream/reshape
%       sym/reshape
%       symbolic/reshape
%       tabular/reshape
%       tall/reshape
%