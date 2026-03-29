% :  Colon.
%    J:K  is the same as [J, J+1, ..., J+m], where m = fix(K-J). In the 
%    case where both J and K are integers, this is simply [J, J+1, ..., K]. 
%    This syntax returns an empty matrix if J > K.
% 
%    J:I:K  is the same as [J, J+I, ..., J+m*I], where m = fix((K-J)/I). 
%    This syntax returns an empty matrix when I == 0, I > 0 and J > K, or 
%    I < 0 and J < K.
% 
%    COLON(J,K) is the same as J:K and COLON(J,I,K) is the same as J:I:K.
% 
%    The colon notation can be used to pick out selected rows, columns
%    and elements of vectors, matrices, and arrays.  A(:) is all the
%    elements of A, regarded as a single column. On the left side of an
%    assignment statement, A(:) fills A, preserving its shape from before.
%    A(:,J) is the J-th column of A.  A(J:K) is [A(J),A(J+1),...,A(K)].
%    A(:,J:K) is [A(:,J),A(:,J+1),...,A(:,K)] and so on.
% 
%    The colon notation can be used with a cell array to produce a comma-
%    separated list.  C{:} is the same as C{1},C{2},...,C{end}.  The comma
%    separated list syntax is valid inside () for function calls, [] for
%    concatenation and function return arguments, and inside {} to produce
%    a cell array.  Expressions such as S(:).name produce the comma separated
%    list S(1).name,S(2).name,...,S(end).name for the structure S.
% 
%    For the use of the colon in the FOR statement, See FOR.
%    For the use of the colon in a comma separated list, See VARARGIN.
%
%    Documentation for colon
%       doc colon
%
%    Other uses of colon
%
%       calendarDuration/colon     duration/colon
%       codistributed/colon        gpuArray/colon
%       codistributor1d/colon      rise_dates.dates/colon
%       codistributor2dbc/colon    sym/colon
%       datetime/colon             symbolic/colon
%       distributed/colon
%