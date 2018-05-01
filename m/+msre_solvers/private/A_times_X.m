function D=A_times_X(A,X,locs)
% H1 line
%
% ::
%
%
% Args:
%
% Returns:
%    :
%
% Note:
%
% Example:
%
%    See also:

% computes A*X assuming the zero columns of A are deleted. the nonzero
% columns are given in the logical vector locs
D=A*X(locs,:); 
end

