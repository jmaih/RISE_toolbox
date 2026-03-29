% ^  Matrix power.
%    Z = X^y is X to the y power if y is a scalar and X is square. If y
%    is an integer greater than one, the power is computed by repeated
%    squaring. For other values of y the calculation involves
%    eigenvalues and eigenvectors.
% 
%    Z = x^Y is x to the Y power if Y is a square matrix and x is a scalar.
%    Computed using eigenvalues and eigenvectors.
% 
%    Z = X^Y, where both X and Y are matrices, is an error.
% 
%    C = MPOWER(A,B) is called for the syntax 'A ^ B' when A or B is an object.
% 
%    See <a href="matlab:helpview('matlab','MATLAB_OPS')">MATLAB Operators and Special Characters</a> for more details.
% 
%    See also POWER.
%
%    Documentation for mpower
%       doc mpower
%
%    Other uses of mpower
%
%       adolm/mpower            gpuArray/mpower     sym/mpower
%       aplanar/mpower          rsymbdiff/mpower    symbolic/mpower
%       codistributed/mpower    splanar/mpower      ts/mpower
%