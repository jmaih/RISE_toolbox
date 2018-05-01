function A=symmetrize(A)
% symmetrize - makes a square matrix symmetric
%
% ::
%
%
%   A=symmetrize(A)
%
% Args:
%
%    - **A** [square matrix]
%
% Returns:
%    :
%
%    - **A** [square matrix]
%
% Note:
%
% Example:
%
%    See also:

A=.5*(A+A.');
end

