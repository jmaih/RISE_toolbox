%  FOLD Apply a binary function left-associatively to a vector
%    FOLD(F, V, INITIALVALUE) iteratively computes 
%    RES = F(V(1), V(2))
%    RES = F(RES, V(3)) 
%    ....
%    RES = F(RES, V(END))
%    and, finally, returns RES. If V is empty, INITIALVALUE is returned. 
%    For non-empty V, the third argument INITIALVALUE may be left out.
% 
%    Example: 
%    syms x
%    fold(@or, x == 1:3) gives x == 1 | x==2 | x==3
%    fold(@and, [], true) gives true   
%    fold(@intersect, {[1, 3], [2, 3, 29], [3, 1, 17]}) gives 3
% 
%    See also PROD, SUM.
%
%    Reference page in Doc Center
%       doc fold
%
%    Other functions named fold
%
%       symfun/fold    ts/fold
%