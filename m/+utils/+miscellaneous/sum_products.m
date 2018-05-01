function res=sum_products(varargin)
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

% sum products of all possible permutations
% example sum_products(A,B,C)=A*B*C + A*C*B + B*A*C + B*C*A + C*A*B + C*B*A
n=length(varargin);
% find all permutations
%----------------------
allperms=cell2mat(utils.gridfuncs.mypermutation(1:n));
res=0;
for irow=1:size(allperms,1)
    % compute product of one permutation
    %-----------------------------------
    tmp=1;
    for icol=1:n
        tmp=tmp*varargin{allperms(irow,icol)};
    end
    % sum all products
    %-----------------
    res=res+tmp;
end
end
