% UNIQUE Set unique.
%    C = UNIQUE(A) for the array A returns the same values as in A but with
%    no repetitions. C will be sorted.
% 
%    C = UNIQUE(A,'rows') for the matrix A returns the unique rows of A.
%    The rows of the matrix C will be in sorted order.
% 
%    [C,IA,IC] = UNIQUE(A) also returns index vectors IA and IC such that
%    C = A(IA) and A = C(IC) (or A(:) = C(IC), if A is a matrix or array).
% 
%    [C,IA,IC] = UNIQUE(A,'rows') also returns index vectors IA and IC such
%    that C = A(IA,:) and A = C(IC,:).
% 
%    [C,IA,IC] = UNIQUE(A,OCCURRENCE) and
%    [C,IA,IC] = UNIQUE(A,'rows',OCCURRENCE) specify which index is returned
%    in IA in the case of repeated values (or rows) in A. The default value
%    is OCCURRENCE = 'first', which returns the index of the first occurrence
%    of each repeated value (or row) in A, while OCCURRENCE = 'last' returns
%    the index of the last occurrence of each repeated value (or row) in A.
% 
%    [C,IA,IC] = UNIQUE(A,'stable') returns the values of C in the same order
%    that they appear in A, while [C,IA,IC] = UNIQUE(A,'sorted') returns the
%    values of C in sorted order. If A is a row vector, then C will be a row
%    vector as well, otherwise C will be a column vector. IA and IC are
%    column vectors. If there are repeated values in A, then IA returns the
%    index of the first occurrence of each repeated value.
% 
%    [C,IA,IC] = UNIQUE(A,'rows','stable') returns the rows of C in the same
%    order that they appear in A, while [C,IA,IC] = UNIQUE(A,'rows','sorted')
%    returns the rows of C in sorted order.
% 
%    The behavior of UNIQUE has changed.  This includes:
%      -	occurrence of indices in IA and IC switched from last to first
%      -	IA and IC will always be column index vectors
% 
%    If this change in behavior has adversely affected your code, you may
%    preserve the previous behavior with:
% 
%         [C,IA,IC] = UNIQUE(A,'legacy')
%         [C,IA,IC] = UNIQUE(A,'rows','legacy')
%         [C,IA,IC] = UNIQUE(A,OCCURRENCE,'legacy')
%         [C,IA,IC] = UNIQUE(A,'rows',OCCURRENCE,'legacy')
% 
%    Examples:
% 
%        a = [9 9 9 9 9 9 8 8 8 8 7 7 7 6 6 6 5 5 4 2 1]
% 
%        [c1,ia1,ic1] = unique(a)
%        % returns
%        c1 = [1 2 4 5 6 7 8 9]
%        ia1 = [21 20 19 17 14 11 7 1]'
%        ic1 = [8 8 8 8 8 8 7 7 7 7 6 6 6 5 5 5 4 4 3 2 1]'
% 
%        [c2,ia2,ic2] = unique(a,'stable')
%        % returns
%        c2 = [9 8 7 6 5 4 2 1]
%        ia2 = [1 7 11 14 17 19 20 21]'
%        ic2 = [1 1 1 1 1 1 2 2 2 2 3 3 3 4 4 4 5 5 6 7 8]'
% 
%        c = unique([1 NaN NaN 2])
%        % NaNs compare as not equal, so this returns
%        c = [1 2 NaN NaN]
% 
%    Class support for input A:
%       - logical, char, all numeric classes
%       - cell arrays of strings
%       -- 'rows' option is not supported for cell arrays
%       - objects with methods SORT (SORTROWS for the 'rows' option) and NE
%       -- including heterogeneous arrays
% 
%    See also UNIQUETOL, INTERSECT, ISMEMBER, ISMEMBERTOL, UNION, SETDIFF,
%             SETXOR, SORT, SORTROWS.
%
%    Documentation for unique
%       doc unique
%
%    Other uses of unique
%
%       calendarDuration/unique    gpuArray/unique
%       categorical/unique         rise_dates.dates/unique
%       cell/unique                sym/unique
%       codistributed/unique       symfun/unique
%       dataset/unique             tabular/unique
%       datetime/unique            tall/unique
%       duration/unique
%