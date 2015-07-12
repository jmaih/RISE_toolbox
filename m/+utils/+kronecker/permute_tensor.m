function varargout=permute_tensor(A1_Ak,matsizes,varargin)
% permute_tensor -- permute a kronecker product
%
% Syntax
% -------
% ::
%
%   [H1,H2,...,Hn]=permute_tensor(A1_Ak,matsizes,order_1,...,order_n)
%
% Inputs
% -------
%
% - **A1_Ak** [matrix]: representing the kronecker product A1*A2*A3*...*Ak
%
% - **matsizes** [matrix]: k x 2 matrix of size of the various matrices
% entering the tensor. Each row represents the size of a matrix and it is
% assumed that the main(or first) tensor product is ordered [1,2,...,k]
%
% - **order_i** [vector]: permutation of [1,2,...,k]. N.B: all elements
% 1,2,...,k should be part of the vector and the vector should have exactly
% k elements
%
% Outputs
% --------
%
% - **Hi** [matrix]: tensor permutation of A1_Ak according to order "i"
%
% More About
% ------------
%
% Examples
% ---------
%
% See also: sum_permutations, tensorperm

persistent A1_Ak0 orig_order nmat rows_prod cols_prod

if isempty(A1_Ak)
    A1_Ak=A1_Ak0;
else
    nmat=size(matsizes,1);
    
    rows_sizes=matsizes(:,1).';
    rows_prod=prod(rows_sizes);
    rowsizflip=fliplr(rows_sizes);
    
    cols_sizes=matsizes(:,2).';
    cols_prod=prod(cols_sizes);
    colsizflip=fliplr(cols_sizes);
    
    orig_order=fliplr(1:nmat);% nmat:-1:1
    
    A1_Ak=reshape(full(A1_Ak),[rowsizflip,colsizflip]);
    A1_Ak0=A1_Ak;
end

varargout=varargin;

for iarg=1:length(varargin)
    new_order=fliplr(varargin{iarg}(:).');
    prows=orig_order(new_order);
    pcols=prows+nmat;
    varargout{iarg}=sparse(...
        reshape(...
        permute(A1_Ak,[prows,pcols]),...
        rows_prod,cols_prod)...
        );
end

end

%{
clc
kronall=@utils.kronecker.kronall;
imax=10;
imin=3;
matrices={
    rand(randi([imin,imax]),randi([imin,imax]))
    rand(randi([imin,imax]),randi([imin,imax]))
    rand(randi([imin,imax]),randi([imin,imax]))
    rand(randi([imin,imax]),randi([imin,imax]))
    };
nmat=numel(matrices);

matsizes=zeros(nmat,2);

for ii=1:nmat
    matsizes(ii,:)=size(matrices{ii});
end

Main=kronall(matrices{:});

orders={randperm(nmat),randperm(nmat),randperm(nmat),randperm(nmat),randperm(nmat)};

no=numel(orders);

Alternatives=cell(1,no);
for ialt=1:no
    tmp=matrices(orders{ialt});
    Alternatives{ialt}=kronall(tmp{:});
end

tic
[new_tensors{1:no}]=utils.kronecker.permute_tensor(Main,matsizes,orders{:});
toc


for ialt=1:no
    max(max(abs(Alternatives{ialt}-new_tensors{ialt})))
end

%}
