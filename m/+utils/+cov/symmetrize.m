function A=symmetrize(A)
% symmetrize - makes a square matrix symmetric
%
% Syntax
% -------
% ::
%
%   A=symmetrize(A)
%
% Inputs
% -------
%
% - **A** [square matrix]
%
% Outputs
% --------
%
% - **A** [square matrix]
%
% More About
% ------------
%
% Examples
% ---------
%
% See also:

A=.5*(A+A.');
end

