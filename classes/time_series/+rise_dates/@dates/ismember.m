% ISMEMBER True for set member.
%    LIA = ISMEMBER(A,B) for arrays A and B returns an array of the same
%    size as A containing true where the elements of A are in B and false
%    otherwise.
% 
%    LIA = ISMEMBER(A,B,'rows') for matrices A and B with the same number
%    of columns, returns a vector containing true where the rows of A are
%    also rows of B and false otherwise.
% 
%    [LIA,LOCB] = ISMEMBER(A,B) also returns an array LOCB containing the
%    lowest absolute index in B for each element in A which is a member of
%    B and 0 if there is no such index.
% 
%    [LIA,LOCB] = ISMEMBER(A,B,'rows') also returns a vector LOCB containing
%    the lowest absolute index in B for each row in A which is a member
%    of B and 0 if there is no such index.
% 
%    The behavior of ISMEMBER has changed.  This includes:
%      -	occurrence of indices in LOCB switched from highest to lowest
%      -	tighter restrictions on combinations of classes
% 
%    If this change in behavior has adversely affected your code, you may
%    preserve the previous behavior with:
% 
%       [LIA,LOCB] = ISMEMBER(A,B,'legacy')
%       [LIA,LOCB] = ISMEMBER(A,B,'rows','legacy')
% 
%    Examples:
% 
%       a = [9 9 8 8 7 7 7 6 6 6 5 5 4 4 2 1 1 1]
%       b = [1 1 1 3 3 3 3 3 4 4 4 4 4 9 9 9]
% 
%       [lia1,locb1] = ismember(a,b)
%       % returns
%       lia1 = [1 1 0 0 0 0 0 0 0 0 0 0 1 1 0 1 1 1]
%       locb1 = [14 14 0 0 0 0 0 0 0 0 0 0 9 9 0 1 1 1]
% 
%       [lia,locb] = ismember([1 NaN 2 3],[3 4 NaN 1])
%       % NaNs compare as not equal, so this returns
%       lia = [1 0 0 1], locb = [4 0 0 1]
% 
%    Class support for inputs A and B, where A and B must be of the same
%    class unless stated otherwise:
%       - logical, char, all numeric classes (may combine with double arrays)
%       - cell arrays of strings (may combine with char arrays)
%       -- 'rows' option is not supported for cell arrays
%       - objects with methods SORT (SORTROWS for the 'rows' option), EQ and NE
%       -- including heterogeneous arrays derived from the same root class
% 
%    See also ISMEMBERTOL, INTERSECT, UNION, UNIQUE, UNIQUETOL, SETDIFF,
%             SETXOR, SORT, SORTROWS.
%
%    Documentation for ismember
%       doc ismember
%
%    Other uses of ismember
%
%       calendarDuration/ismember    gpuArray/ismember
%       categorical/ismember         mtree/ismember
%       cell/ismember                rise_dates.dates/ismember
%       codistributed/ismember       sym/ismember
%       dataset/ismember             tabular/ismember
%       datetime/ismember            tall/ismember
%       duration/ismember
%