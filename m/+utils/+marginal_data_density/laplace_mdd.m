function log_mdd=laplace_mdd(log_post,Hinv)

% LAPLACE_MDD -- computes the log marginal data density using the laplace
% approximation
%
% Syntax
% -------
% ::
%
%   log_mdd=LAPLACE_MDD(log_post,Hinv)
%
% Inputs
% -------
%
% - **log_post** [numeric]: log of posterior kernel evaluated at the mode
%
% - **Hinv** [matrix]: inverse Hessian (negative of the second derivatives
% of the log posterior kernel)
%
% Outputs
% --------
%
% - **log_mdd** [numeric]: log marginal data density
%
% More About
% ------------
%
% Examples
% ---------
%
% See also: 

d=size(Hinv,1);

% log_mdd=.5*npar*log(2*pi)-.5*log(det(H))+log_post;
% or alternatively

log_mdd=.5*d*log(2*pi)+.5*log(det(Hinv))+log_post;

end