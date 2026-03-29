% INTERSECT Set intersection.
%    C = INTERSECT(A,B) for vectors A and B, returns the values common to
%    the two vectors with no repetitions. C will be sorted.
% 
%    C = INTERSECT(A,B,'rows') for matrices A and B with the same
%    number of columns, returns the rows common to the two matrices. The
%    rows of the matrix C will be in sorted order.
% 
%    [C,IA,IB] = INTERSECT(A,B) also returns index vectors IA and IB such
%    that C = A(IA) and C = B(IB). If there are repeated common values in
%    A or B then the index of the first occurrence of each repeated value is
%    returned.
% 
%    [C,IA,IB] = INTERSECT(A,B,'rows') also returns index vectors IA and IB
%    such that C = A(IA,:) and C = B(IB,:).
% 
%    [C,IA,IB] = INTERSECT(A,B,'stable') for arrays A and B, returns the
%    values of C in the same order that they appear in A.
%    [C,IA,IB] = INTERSECT(A,B,'sorted') returns the values of C in sorted
%    order.
%    If A and B are row vectors, then C will be a row vector as well,
%    otherwise C will be a column vector. IA and IB are column vectors.
%    If there are repeated common values in A or B then the index of the
%    first occurrence of each repeated value is returned.
% 
%    [C,IA,IB] = INTERSECT(A,B,'rows','stable') returns the rows of C in the
%    same order that they appear in A.
%    [C,IA,IB] = INTERSECT(A,B,'rows','sorted') returns the rows of C in
%    sorted order.
% 
%    The behavior of INTERSECT has changed.  This includes:
%      -	occurrence of indices in IA and IB switched from last to first
%      -	orientation of vector C
%      -	IA and IB will always be column index vectors
%      -	tighter restrictions on combinations of classes
% 
%    If this change in behavior has adversely affected your code, you may
%    preserve the previous behavior with:
% 
%       [C,IA,IB] = INTERSECT(A,B,'legacy')
%       [C,IA,IB] = INTERSECT(A,B,'rows','legacy')
% 
%    Examples:
% 
%       a = [9 9 9 9 9 9 8 8 8 8 7 7 7 6 6 6 5 5 4 2 1]
%       b = [1 1 1 3 3 3 3 3 4 4 4 4 4 10 10 10]
% 
%       [c1,ia1,ib1] = intersect(a,b)
%       % returns
%       c1 = [1 4], ia1 = [21 19]', ib1 = [1 9]'
% 
%       [c2,ia2,ib2] = intersect(a,b,'stable')
%       % returns
%       c2 = [4 1], ia2 = [19 21]', ib2 = [9 1]'
% 
%       c = intersect([1 NaN 2 3],[3 4 NaN 1])
%       % NaNs compare as not equal, so this returns
%       c = [1 3]
% 
%    Class support for inputs A and B, where A and B must be of the same
%    class unless stated otherwise:
%       - logical, char, all numeric classes (may combine with double arrays)
%       - cell arrays of strings (may combine with char arrays)
%       -- 'rows' option is not supported for cell arrays
%       - objects with methods SORT (SORTROWS for the 'rows' option), EQ and NE
%       -- including heterogeneous arrays derived from the same root class
% 
%    See also UNIQUE, UNION, SETDIFF, SETXOR, ISMEMBER, SORT, SORTROWS.
%
%    Documentation for intersect
%       doc intersect
%
%    Other uses of intersect
%
%       calendarDuration/intersect    gpuArray/intersect
%       categorical/intersect         polyshape/intersect
%       cell/intersect                rise_dates.dates/intersect
%       codistributed/intersect       tabular/intersect
%       dataset/intersect             tall/intersect
%       datetime/intersect            ts/intersect
%       duration/intersect
%