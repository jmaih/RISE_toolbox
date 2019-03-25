% SUBSASGN Subscripted assignment.
%    A(I) = B assigns the values of B into the elements of A specified by
%    the subscript vector I.  B must have the same number of elements as I
%    or be a scalar. 
% 
%    A(I,J) = B assigns the values of B into the elements of the rectangular
%    submatrix of A specified by the subscript vectors I and J.  A colon used as
%    a subscript, as in A(I,:) = B, indicates all columns of those rows
%    indicated by vector I. Similarly, A(:,J) = B means all rows of columns J.
% 
%    A(I,J,K,...) = B assigns the values of B to the submatrix of A specified
%    by the subscript vectors I, J, K, etc. A colon used as a subscript, as in
%    A(I,:,K) = B, indicates the entire dimension. 
% 
%    For both A(I,J) = B and the more general multi-dimensional 
%    A(I,J,K,...) = B, B must be LENGTH(I)-by-LENGTH(J)-by-LENGTH(K)-... , or
%    be shiftable to that size by adding or removing singleton dimensions, or
%    contain a scalar, in which case its value is replicated to form a matrix
%    of that size.
% 
%    A{I} = B when A is a cell array and I is a scalar places a copy of
%    the array B into the specified cell of A.  If I has more than one
%    element, this expression is an error.  Use [A{I}] = DEAL(B) to place
%    copies of B into multiple cells of A.  Multiple subscripts that
%    specify a scalar element, as in A{3,4} = magic(3), also work.
% 
%    A(I).field = B when A is a structure array and I is a scalar places
%    a copy of the array B into the field with the name 'field'.  If I
%    has more than one element, this expression is an error. Use
%    [A(I).field] = DEAL(B) to place copies of B into multiple elements
%    of A.  If A is a 1-by-1 structure array, then the subscript can be
%    dropped.  In this case, A.field = B is the same as A(1).field = B.
% 
%    A = SUBSASGN(A,S,B) is called for the syntax A(I)=B, A{I}=B, or
%    A.I=B when A is an object.  S is a structure array with the fields:
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
%    S.type='{}' and the syntax A.field=B invokes SUBSASGN(A,S,B) where
%    S.type='.' and S.subs='field'.
% 
%    These simple calls are combined in a straightforward way for
%    more complicated subscripting expressions.  In such cases
%    length(S) is the number of subscripting levels.  For instance, 
%    A(1,2).name(3:5)=B invokes A=SUBSASGN(A,S,B) where S is 3-by-1
%    structure array with the following values:
%        S(1).type='()'       S(2).type='.'        S(3).type='()'
%        S(1).subs={1,2}      S(2).subs='name'     S(3).subs={3:5}
% 
%    See also SUBSREF, SUBSTRUCT, PAREN, SUBSINDEX, LISTS, 
%             NUMARGUMENTSFROMSUBSCRIPT.
%
%    Reference page in Doc Center
%       doc subsasgn
%
%    Other functions named subsasgn
%
%       axisobj/subsasgn             gbm/subsasgn
%       axistext/subsasgn            gpuArray/subsasgn
%       bm/subsasgn                  graph/subsasgn
%       calendarDuration/subsasgn    heston/subsasgn
%       categorical/subsasgn         hgbin/subsasgn
%       cev/subsasgn                 hwv/subsasgn
%       cir/subsasgn                 instrument/subsasgn
%       classregtree/subsasgn        scribehandle/subsasgn
%       codistributed/subsasgn       scribehgobj/subsasgn
%       Composite/subsasgn           sde/subsasgn
%       dataset/subsasgn             sdeddo/subsasgn
%       datetime/subsasgn            sdeld/subsasgn
%       diffusion/subsasgn           sdemrd/subsasgn
%       digraph/subsasgn             serial/subsasgn
%       distributed/subsasgn         sym/subsasgn
%       drift/subsasgn               tabular/subsasgn
%       duration/subsasgn            tall/subsasgn
%       fighandle/subsasgn           ts/subsasgn
%       figobj/subsasgn              tscollection/subsasgn
%       fints/subsasgn
%