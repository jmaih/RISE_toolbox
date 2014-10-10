function res=sum_kroneckers(varargin)
% H1 line
%
% Syntax
% -------
% ::
%
% Inputs
% -------
%
% Outputs
% --------
%
% Description
% ------------
%
% Examples
% ---------
%
% See also: 

% sum products of all possible permutations
% example sum_kroneckers(A,B,C)=kron(A,B,C) + kron(A,C,B) + kron(B,A,C) +
% kron(B,C,A) + kron(C,A,B) + kron(C,B,A)
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
