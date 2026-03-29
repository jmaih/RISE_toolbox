% SUBSASGN Subscripted assignment.
%    A(I) = B assigns the values of B into the elements of A specified by
%    the subscript vector I.  B must have the same number of elements as I
%    or be a scalar. 
% 
%    A(I,J) = B assigns the values of B into the elements of the rectangular
%    submatrix of A specified by the subscript vectors I and J.  A colon
%    used as a subscript, as in A(I,:) = B, indicates all columns of those
%    rows indicated by vector I. Similarly, A(:,J) = B means all rows of
%    columns J.
% 
%    A(I,J,K,...) = B assigns the values of B to the subarray of A specified
%    by the subscript vectors I, J, K, etc. A colon used as a subscript, as
%    in A(I,:,K) = B, indicates the entire dimension.
% 
%    For both A(I,J) = B and the more general multi-dimensional A(I,J,K,...)
%    = B, B must be LENGTH(I)-by-LENGTH(J)-by-LENGTH(K)-... , or be
%    shiftable to that size by adding or removing singleton dimensions, or
%    contain a scalar, in which case its value is replicated to form a
%    matrix of that size.
% 
%    A{I} = B when A is a cell array and I is a scalar places a copy of the
%    array B into the specified cell of A.  If I has more than one element,
%    this expression is an error.  Use [A{I}] = DEAL(B) to place copies of B
%    into multiple cells of A.  Multiple subscripts that specify a scalar
%    element, as in A{3,4} = magic(3), also work.
% 
%    A(I).label = B when A is a structure or object array and I is a scalar
%    places a copy of the array B into the field or property with the name
%    'label'.  If I has more than one element, this expression results in an
%    error. Use [A(I).label] = DEAL(B) to place copies of B into fields or
%    properties of multiple elements of A. If A is a 1-by-1 array, then the
%    subscript can be dropped.  In this case, A.label = B is the same as
%    A(1).label = B.
% 
%    When var is a variable containing 'label', A(I).(var) = B places a copy
%    of the array B into the field or property with the name 'label'. When I
%    has more than one element, use [A(I).(var)] = DEAL(B) to place copies
%    of B into fields or properties of multiple elements of A.
% 
%    A class can implement a method called SUBSASGN to overload indexed
%    assignment. When A is an instance of a class that has a SUBSASGN
%    method, A = SUBSASGN(A,S,B) is called for the syntax A(I)=B, A{I}=B,
%    A.I=B, or A.(I)=B. The argument S is a structure array with the fields:
%        type -- character vector or string containing '()', '{}', or '.' 
%                specifying the subscript type.
%        subs -- cell array, character vector, or string containing the 
%                actual subscripts.
% 
%    For instance, the syntax A(1:2,:)=B calls A=SUBSASGN(A,S,B) where
%    S is a 1-by-1 structure with S.type='()' and S.subs = {1:2,':'}. A
%    colon used as a subscript is passed as the character vector ':'.
% 
%    Similarly, the syntax A{1:2}=B invokes A=SUBSASGN(A,S,B) where
%    S.type='{}' and S.subs={1:2}. The syntax A.field=B invokes
%    SUBSASGN(A,S,B) where S.type='.' and S.subs='field'.
% 
%    These simple calls are combined in a straightforward way for more
%    complicated subscripting expressions.  In such cases length(S) is the
%    number of subscripting levels.  For instance, A(1,2).name(3:5)=B
%    invokes A=SUBSASGN(A,S,B) where S is a 1-by-3 structure array with the
%    following values:
%        S(1).type='()'       S(2).type='.'        S(3).type='()'
%        S(1).subs={1,2}      S(2).subs='name'     S(3).subs={3:5}
% 
%    Instead of implementing SUBSASGN, class authors can use the mixins in
%    the matlab.mixin.indexing package to overload indexing. See the
%    documentation for more information.
% 
%    See also SUBSREF, SUBSTRUCT, numArgumentsFromSubscript, SUBSINDEX, 
%             PAREN, LISTS, matlab.mixin.indexing
%
%    Documentation for subsasgn
%       doc subsasgn
%
%    Other uses of subsasgn
%
%       axisobj/subsasgn             graph/subsasgn
%       axistext/subsasgn            heston/subsasgn
%       bates/subsasgn               hgbin/subsasgn
%       bm/subsasgn                  hwv/subsasgn
%       calendarDuration/subsasgn    instrument/subsasgn
%       categorical/subsasgn         merton/subsasgn
%       cev/subsasgn                 parallel-computing/subsasgn
%       cir/subsasgn                 RandStream/subsasgn
%       codistributed/subsasgn       rvm/subsasgn
%       Composite/subsasgn           scribehandle/subsasgn
%       dataset/subsasgn             scribehgobj/subsasgn
%       datetime/subsasgn            sde/subsasgn
%       diffusion/subsasgn           sdeddo/subsasgn
%       digraph/subsasgn             sdeld/subsasgn
%       distributed/subsasgn         sdemrd/subsasgn
%       drift/subsasgn               serial/subsasgn
%       duration/subsasgn            sym/subsasgn
%       fighandle/subsasgn           tall/subsasgn
%       figobj/subsasgn              ts/subsasgn
%       gbm/subsasgn                 tscollection/subsasgn
%