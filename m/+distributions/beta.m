% BETA   Beta function.
%    Y = BETA(Z,W) computes the beta function for corresponding
%    elements of Z and W.  The beta function is defined as
% 
%    beta(z,w) = integral from 0 to 1 of t.^(z-1) .* (1-t).^(w-1) dt.
% 
%    The arrays Z and W must be real and nonnegative. Both arrays must be 
%    the same size, or either can be scalar.
% 
%    Class support for inputs Z,W:
%       float: double, single
% 
%    See also BETAINC, BETALN.
%
%    Reference page in Doc Center
%       doc beta
%
%    Other functions named beta
%
%       codistributed/beta    gpuArray/beta    sym/beta
%