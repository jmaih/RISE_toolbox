function D=X_times_A(X,A,locs)
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

% computes X*A assuming the zero columns of A are deleted. the nonzero
% columns are given in the logical vector locs
D=zeros(size(X,1),size(A,1));
D(:,locs)=X*A; % <--- D(:,locs)=X*A(:,locs);
end
