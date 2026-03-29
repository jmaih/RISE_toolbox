%--- parse.m not found. Showing help for sparse instead. ---
%
% sparse - Create sparse matrix
%    This MATLAB function converts a full matrix into sparse form by
%    squeezing out any zero elements.
%
%    Syntax
%      S = sparse(A)
%
%      S = sparse(m,n)
%
%      S = sparse(i,j,v)
%      S = sparse(i,j,v,m,n)
%      S = sparse(i,j,v,m,n,nz)
%
%    Input Arguments
%      A - Input matrix
%        full matrix | sparse matrix
%      i,j - Subscript pairs (as separate arguments)
%        scalars | vectors | matrices
%      v - Values
%        scalar | vector | matrix
%      m,n - Size of each dimension (as separate arguments)
%        integer values
%      nz - Storage allocation for nonzero elements
%        nonnegative integer
%
%    Examples
%      openExample('matlab/ConvertFullMatrixToSparseStorageExample')
%      openExample('matlab/CreateAllZeroSparseMatrixExample')
%      openExample('matlab/CreateSparseMatrixOfNonzerosWithSpecifiedSizeExample')
%      openExample('matlab/PreallocateStorageForNonzerosInSparseMatrixExample')
%      openExample('matlab/AccumulateValuesIntoSparseMatrixExample')
%
%    See also diag, find, full, issparse, nnz, nonzeros, nzmax, spones,
%      sprandn, sprandsym, spy, accumarray, spalloc, speye
%
%    Introduced in MATLAB before R2006a
%    Documentation for sparse
%       doc sparse
%
%    Other uses of sparse
%
%       codistributed/sparse      codistributor2dbc/sparse
%       codistributor1d/sparse    gpuArray/sparse
%