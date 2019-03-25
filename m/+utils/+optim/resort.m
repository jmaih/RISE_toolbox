%--- help for sort ---
%
% SORT   Sort in ascending or descending order.
%    B = SORT(A) sorts in ascending order.
%    The sorted output B has the same type and size as A:
%    - For vectors, SORT(A) sorts the elements of A in ascending order.
%    - For matrices, SORT(A) sorts each column of A in ascending order.
%    - For N-D arrays, SORT(A) sorts along the first non-singleton dimension.
% 
%    B = SORT(A,DIM) also specifies a dimension DIM to sort along.
% 
%    B = SORT(A,DIRECTION) and B = SORT(A,DIM,DIRECTION) also specify the
%    sort direction. DIRECTION must be:
%        'ascend'  - (default) Sorts in ascending order.
%        'descend' - Sorts in descending order.
% 
%    B = SORT(A,...,'MissingPlacement',M) also specifies where to place the
%    missing elements (NaN/NaT/<undefined>/<missing>) of A. M must be:
%        'auto'  - (default) Places missing elements last for ascending sort
%                  and first for descending sort.
%        'first' - Places missing elements first.
%        'last'  - Places missing elements last.
% 
%    B = SORT(A,...,'ComparisonMethod',C) specifies how to sort complex
%    numbers. The comparison method C must be:
%        'auto' - (default) Sorts real numbers according to 'real', and
%                 complex numbers according to 'abs'.
%        'real' - Sorts according to REAL(A). Elements with equal real parts
%                 are then sorted by IMAG(A).
%        'abs'  - Sorts according to ABS(A). Elements with equal magnitudes
%                 are then sorted by ANGLE(A).
% 
%    [B,I] = SORT(A,...) also returns a sort index I which specifies how the
%    elements of A were rearranged to obtain the sorted output B:
%    - If A is a vector, then B = A(I).  
%    - If A is an m-by-n matrix and DIM = 1, then
%        for j = 1:n, B(:,j) = A(I(:,j),j); end
% 
%    The sort ordering is stable. Namely, when more than one element has the
%    same value, the order of the equal elements is preserved in the sorted
%    output B and the indices I relating to equal elements are ascending.
% 
%    Examples:
%      % Sort a vector in ascending order
%        sort([0 3 1 0 2 0 1 6])
%      % Sort each column or row of a matrix
%        A = [3 7 5; 0 4 2]
%        B1 = sort(A,1)  % sort each column
%        B2 = sort(A,2)  % sort each row
%      % Sort complex numbers according to their real part
%        A = [1+1j ;  1-1j ;  -2-2j ;  0 ;  -2+2j]
%        B = sort(A,'ComparisonMethod','real')
% 
%    See also ISSORTED, SORTROWS, MIN, MAX, MINK, MAXK.
%
%    Reference page in Doc Center
%       doc sort
%
%