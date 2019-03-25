%--- help for vectorize ---
%
% VECTORIZE Vectorize expression.
%    VECTORIZE(S), when S is a string expression, inserts a '.' before
%    any '^', '*' or '/' in S.  The result is a character string.
% 
%    VECTORIZE will not accept INLINE function objects in a future
%    release. Use anonymous functions and FUNC2STR instead.
% 
%    VECTORIZE(FUN), when FUN is an INLINE function object, vectorizes the
%    formula for FUN.  The result is the vectorized version of the INLINE
%    function.
% 
%    See also INLINE/FORMULA, INLINE, FUNCTION_HANDLE.
%
%    Reference page in Doc Center
%       doc vectorize
%
%