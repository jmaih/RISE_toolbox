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
%    A{I} when A is a cell array and I is a scalar is a copy of the array in
%    the specified cell of A.  If I has more than one element, this
%    expression is a comma separated list (see LISTS). Multiple subscripts
%    that specify a scalar element, as in A{3,4}, also work.
% 
%    A(I).label when A is a structure or object array and I is a scalar is a
%    copy of the array in the field or property with the name 'label'. If I
%    has more than one element, this expression is a comma separated list.
%    If A is a 1-by-1 array, then the subscript can be dropped. In this
%    case, A.label is the same as A(1).label.
% 
%    When var is a variable containing 'label', A(I).(var) is a copy of the
%    array in the field or property with the name 'label'. If I has more
%    than one element, this expression is a comma separated list.
% 
%    A class can implement a method named SUBSREF to overload indexed
%    reference. When a class has a SUBSREF method, B = SUBSREF(A,S) is
%    called for the syntax A(I), A{I}, A.I, or A.(I) whenever A is an
%    instance of the class. The argument S is a structure array with the
%    fields:
%        type -- character vector or string containing '()', '{}', or '.' 
%                specifying the subscript type.
%        subs -- cell array, character vector, or string containing the 
%                actual subscripts.
% 
%    Subscripting expressions can use more than one level to form more
%    complicated expressions, for example A{1}.field(3:5). For an expression
%    with N subscripting levels, S is a 1-by-N structure array.
%    
%    Instead of implementing SUBSREF, class authors can use the mixins in
%    the matlab.mixin.indexing package to overload indexing. See the
%    documentation for more information.
% 
%    See also SUBSASGN, SUBSTRUCT, numArgumentsFromSubscript, SUBSINDEX,
%             PAREN, LISTS, matlab.mixin.indexing
%
%    Documentation for subsref
%       doc subsref
%
%    Other uses of subsref
%
%       axisobj/subsref             graph/subsref
%       axistext/subsref            hgbin/subsref
%       calendarDuration/subsref    inline/subsref
%       categorical/subsref         instrument/subsref
%       codistributed/subsref       parallel-computing/subsref
%       Composite/subsref           printtemplate/subsref
%       dataset/subsref             RandStream/subsref
%       datetime/subsref            scribehandle/subsref
%       diffusion/subsref           scribehgobj/subsref
%       digraph/subsref             serial/subsref
%       distributed/subsref         sym/subsref
%       drift/subsref               tall/subsref
%       duration/subsref            ts/subsref
%       fighandle/subsref           tscollection/subsref
%       figobj/subsref
%