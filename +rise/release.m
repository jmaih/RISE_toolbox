%--- help for sym/release ---
%
% RELEASE   Release held integrals.
%    RELEASE(S) removes the 'Hold'-flag from all integrals in S.
% 
%    Example:
%      syms x
%      J = int(x, x, 'Hold', true) returns the symbolic call to int
%      release(J) returns x^2/2
% 
%    See also INT.
%
%    Other uses of release
%
%       COM/release    symbolic/release    system/release
%