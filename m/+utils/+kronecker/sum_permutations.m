function P=sum_permutations(ABCD,matsizes,varargin)
% sum_permutations -- sums the permutations of a tensor product
%
% Syntax
% -------
% ::
%
%   P=sum_permutations(ABCD,matsizes,order_1,order_2,...,order_n)
%
%   P=sum_permutations(ABCD,matsizes,{order_1,order_2,...,order_n})
%
%   P=sum_permutations(ABCD,matsizes,order_1,order_2,...,order_n,'grid')
%
%   P=sum_permutations(ABCD,matsizes,{order_1,order_2,...,order_n},'grid')
%
%   P=sum_permutations(ABCD,matsizes,skip_first,order_1,order_2,...,order_n)
%
%   P=sum_permutations(ABCD,matsizes,skip_first,{order_1,order_2,...,order_n})
%
%   P=sum_permutations(ABCD,matsizes,skip_first,order_1,order_2,...,order_n,'grid')
%
%   P=sum_permutations(ABCD,matsizes,skip_first,{order_1,order_2,...,order_n},'grid')
%
% Inputs
% -------
%
% - **ABCD** [matrix]: tensor product of k matrices
%
% - **matsizes** [matrix]: k x 2 matrix of size of the various matrices
% entering the tensor. Each row represents the size of a matrix and it is
% assumed that the main(or first) tensor product is ordered [1,2,...,k]
%
% - **skip_first** [true|{false}]: if true, the original input matrix is
% not added to the sum
%
% - **order_i** [vector]: permutation of [1,2,...,k]. N.B: all elements
% 1,2,...,k should be part of the vector and the vector should have exactly
% k elements
%
% - **grid** [string]: if present, a grid is used to compute the indexes of
% the main kronecker product
%
% Outputs
% --------
%
% - **P** sum of ABCD and its permutations
%
% More About
% ------------
%
% Examples
% ---------
%
% See also: tensorperm

skip_first=false;
if islogical(varargin{1})
    skip_first=varargin{1};
    varargin=varargin(2:end);
end

is_grid = ischar(varargin{end});

if is_grid
    varargin=varargin(1:end-1);
end

if length(varargin)==1 && iscell(varargin{1})
    orders=varargin{1};
else
    orders=varargin;
end

P=ABCD;
if skip_first
    P=0*P;
end

if is_grid
    [irows,jcols]=utils.kronecker.tensorperm(matsizes,orders{:},'grid');
else
    [irows,jcols]=utils.kronecker.tensorperm(matsizes,orders{:});
end

for ii=1:size(irows,2)
    P=P+ABCD(irows(:,ii),jcols(:,ii));
end

end
