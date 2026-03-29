% HORZCAT Horizontal concatenation.
%    [A B] is the horizontal concatenation of matrices A and B.  A and B
%    must have the same number of rows.  [A,B] is the same thing.  Any
%    number of matrices can be concatenated within one pair of brackets.
%    Horizontal and vertical concatenation can be combined together as in
%    [1 2;3 4].  
% 
%    [A B; C] is allowed if the number of rows of A equals the number of
%    rows of B and the number of columns of A plus the number of columns
%    of B equals the number of columns of C.  The matrices in a
%    concatenation expression can themselves by formed via a
%    concatenation as in [A B;[C D]].  These rules generalize in a
%    hopefully obvious way to allow fairly complicated constructions.
% 
%    N-D arrays are concatenated along the second dimension. The first and
%    remaining dimensions must match.
% 
%    C = HORZCAT(A,B) is called for the syntax '[A  B]' when A or B is an
%    object.
% 
%    Y = HORZCAT(X1,X2,X3,...) is called for the syntax '[X1 X2 X3 ...]'
%    when any of X1, X2, X3, etc. is an object.
% 
%    See also VERTCAT, CAT.
%
%    Documentation for horzcat
%       doc horzcat
%
%    Other uses of horzcat
%
%       calendarDuration/horzcat
%       categorical/horzcat
%       codistributed/horzcat
%       dataset/horzcat
%       datetime/horzcat
%       double/horzcat
%       duration/horzcat
%       gpuArray/horzcat
%       inline/horzcat
%       matlab.mixin.Heterogeneous/horzcat
%       matlab.mixin.indexing.RedefinesParen/horzcat
%       RandStream/horzcat
%       serial/horzcat
%       sym/horzcat
%       tabular/horzcat
%       tall/horzcat
%       ts/horzcat
%       tscollection/horzcat
%