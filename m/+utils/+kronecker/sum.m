function res=sum(varargin)
% sum sums products of all possible permutations of kronecker products
%
% ::
%
%
%   res=sum(A1,A2,...,Ak)
%
% Args:
%
%    - **Ai** [matrix]: matrices for which we want to take the sum of the
%      kronecker product
%
% Returns:
%    :
%
%    - **res** [matrix]: sum of the kronecker products
%
% Note:
%
%    sum(A,B,C)=kron(A,B,C) + kron(A,C,B) + kron(B,A,C) + kron(B,C,A) +
%    kron(C,A,B) + kron(C,B,A)
%
% Example:
%
%    See also:

n=length(varargin);
% find all permutations
%----------------------
allperms=cell2mat(utils.gridfuncs.mypermutation(1:n));
res=0;
for irow=1:size(allperms,1)
    % compute kroneckers of one permutation
    %--------------------------------------
    tmp=varargin{allperms(irow,1)};
    for icol=2:n
        tmp=kron(tmp,varargin{allperms(irow,icol)});
    end
    % sum all products
    %-----------------
    res=res+tmp;
end
end
