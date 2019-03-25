% /   Slash or right matrix divide.
%    A/B is the matrix division of B into A, which is roughly the
%    same as A*INV(B) , except it is computed in a different way.
%    More precisely, A/B = (B'\A')'. See MLDIVIDE for details.
% 
%    C = MRDIVIDE(A,B) is called for the syntax 'A / B' when A or B is an
%    object.
% 
%    See also MLDIVIDE, RDIVIDE, LDIVIDE.
%
%    Reference page in Doc Center
%       doc mrdivide
%
%    Other functions named mrdivide
%
%       codistributed/mrdivide    gpuArray/mrdivide    tall/mrdivide
%       distributed/mrdivide      LagOp/mrdivide       timeseries/mrdivide
%       duration/mrdivide         sym/mrdivide         ts/mrdivide
%       fints/mrdivide
%