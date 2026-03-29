% / Right matrix divide.
%    B/A is the matrix division of A into B, which is roughly the
%    same as B*INV(A) , except it is computed in a different way.
%    More precisely, B/A = (A'\B')'. See MLDIVIDE for details.
% 
%    C = MRDIVIDE(B,A) is called for the syntax 'B / A' when B or A is an object.
% 
%    See <a href="matlab:helpview('matlab','MATLAB_OPS')">MATLAB Operators and Special Characters</a> for more details.
% 
%    See also MLDIVIDE, RDIVIDE, LDIVIDE, PAGEMRDIVIDE.
%
%    Documentation for mrdivide
%       doc mrdivide
%
%    Other uses of mrdivide
%
%       adolm/mrdivide            duration/mrdivide     sym/mrdivide
%       aplanar/mrdivide          gpuArray/mrdivide     symbolic/mrdivide
%       codistributed/mrdivide    LagOp/mrdivide        tall/mrdivide
%       decomposition/mrdivide    rsymbdiff/mrdivide    timeseries/mrdivide
%       distributed/mrdivide      splanar/mrdivide      ts/mrdivide
%