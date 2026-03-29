% PERMUTE Permute array dimensions.
%    B = PERMUTE(A,ORDER) rearranges the dimensions of A so that they are in
%    the order specified by the vector ORDER.  The resulting array has the
%    same values as A but the order of the subscripts needed to access any
%    particular element is rearranged as specified by ORDER.  For an N-D
%    array A, numel(ORDER)>=ndims(A).  All the elements of ORDER must be
%    unique.
%  
%    PERMUTE and IPERMUTE are a generalization of transpose (.') 
%    for N-D arrays.
% 
%    Example:
%       a = rand(1,2,3,4);
%       size(permute(a,[3 2 1 4])) % now it's 3-by-2-by-1-by-4.
% 
%    See also IPERMUTE, CIRCSHIFT, SIZE, PAGETRANSPOSE.
%
%    Documentation for permute
%       doc permute
%
%    Other uses of permute
%
%       calendarDuration/permute    duration/permute      tabular/permute
%       categorical/permute         gpuArray/permute      tall/permute
%       codistributed/permute       RandStream/permute    ts/permute
%       datetime/permute            sym/permute
%