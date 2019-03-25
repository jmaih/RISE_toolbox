% SUBSREF Subscripted reference.
%    A(I) is an array formed from the elements of A specified by the
%    subscript vector I.  The resulting array is the same size as I except
%    for the special case where A and I are both vectors.  In this case,
%    A(I) has the same number of elements as I but has the orientation of A.
% 
%    A(I,J) is an array formed from the elements of the rectangular
%    submatrix of A specified by the subscript vectors I and J.  The
%    resulting array has LENGTH(I) rows and LENGTH(J) columns.  A colon used
%    as a subscript, as in A(I,:), indicates all columns of those rows
%    indicated by vector I. Similarly, A(:,J) = B means all rows of columns
%    J.
% 
%    For multi-dimensional arrays, A(I,J,K,...) is the subarray specified by
%    the subscripts.  The result is LENGTH(I)-by-LENGTH(J)-by-LENGTH(K)-...
% 
%    A{I} when A is a cell array and I is a scalar is a copy of
%    the array in the specified cell of A.  If I has more than one
%    element, this expression is a comma separated list (see LISTS).
%    Multiple subscripts that specify a scalar element, as in A{3,4}, also
%    work.
% 
%    A(I).field when A is a structure array and I is a scalar is a copy of
%    the array in the field with the name 'field'.  If I has more than one
%    element, this expression is a comma separated list.  If A is a 1-by-1
%    structure array, then the subscript can be dropped.  In this case,
%    A.field is the same as A(1).field.
% 
%    B = SUBSREF(A,S) is called for the syntax A(I), A{I}, or A.I when A is
%    an object.  S is a structure array with the fields:
%        type -- character vector or string containing '()', '{}', or '.' 
%                specifying the subscript type.
%        subs -- cell array, character vector, or string containing the 
%                actual subscripts.
% 
%    Subscripting expressions can use more than one level to form more
%    complicated expressions, for example A{1}.field(3:5).  For an
%    expression with N subscripting levels, S is an N-by-1 structure array.
%    See the subscripting examples in the documentation for more
%    information.
% 
%    See also SUBSASGN, SUBSTRUCT, PAREN, SUBSINDEX, LISTS, 
%             NUMARGUMENTSFROMSUBSCRIPT.
%
%    Reference page in Doc Center
%       doc subsref
%
%    Other functions named subsref
%
%       axisobj/subsref             fints/subsref
%       axistext/subsref            gpuArray/subsref
%       calendarDuration/subsref    graph/subsref
%       categorical/subsref         hgbin/subsref
%       classregtree/subsref        inline/subsref
%       codistributed/subsref       instrument/subsref
%       Composite/subsref           printtemplate/subsref
%       dataset/subsref             scribehandle/subsref
%       datetime/subsref            scribehgobj/subsref
%       diffusion/subsref           serial/subsref
%       digraph/subsref             sym/subsref
%       distributed/subsref         tabular/subsref
%       drift/subsref               tall/subsref
%       duration/subsref            ts/subsref
%       fighandle/subsref           tscollection/subsref
%       figobj/subsref
%