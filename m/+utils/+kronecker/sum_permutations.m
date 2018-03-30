function P=sum_permutations(ABCD,matsizes,options,varargin)
% sum_permutations -- sums the permutations of a tensor product
%
% Syntax
% -------
% ::
%
%   P=sum_permutations(ABCD,matsizes,[],order_1,order_2,...,order_n)
%
%   P=sum_permutations(ABCD,matsizes,[],{order_1,order_2,...,order_n})
%
%   P=sum_permutations(ABCD,matsizes,options,order_1,order_2,...,order_n)
%
%   P=sum_permutations(ABCD,matsizes,options,{order_1,order_2,...,order_n})
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
% - **options** [empty|struct]: structure with various options such as
%   - **use_old_algo** [true|{false}]: old and potentially slow algorithm
%   - **use_grid** [true|{false}]: use grid in the old algorithm:a grid is
%       used to compute the indexes of the main kronecker product 
%   - **skip_first** [true|{false}]: if true, the original input matrix is
%   not added to the sum
%
% - **order_i** [vector]: permutation of [1,2,...,k]. N.B: all elements
% 1,2,...,k should be part of the vector and the vector should have exactly
% k elements
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

default_options={
    'use_old_algo',false,@(x)islogical(x),'use_old_algo must be a logical'
    'use_grid',false,@(x)islogical(x),'use_old_algo must be a logical'
    'skip_first',false,@(x)islogical(x),'use_old_algo must be a logical'
    };

if isempty(options)
    
    options=cell2struct(default_options(:,2),default_options(:,1),1);
    
else
    
    if ~isstruct(options)
        
        error('options must be a structure or empty')
        
    end
    
    options=parse_arguments(default_options,options);
    
end

if length(varargin)==1 && iscell(varargin{1})
    
    orders=varargin{1};
    
else
    
    orders=varargin;
    
end

P=ABCD;

if options.skip_first
    
    P=0*P;
    
end

if options.use_old_algo
    
    if options.use_grid
        
        [irows,jcols]=utils.kronecker.tensorperm(matsizes,orders{:},'grid');
        
    else
        
        [irows,jcols]=utils.kronecker.tensorperm(matsizes,orders{:});
        
    end
    
    for ii=1:size(irows,2)
        
        P=P+ABCD(irows(:,ii),jcols(:,ii));
        
    end
    
else
    
    for ii=1:numel(orders)
        
        P=P+utils.kronecker.permute_tensor(ABCD,matsizes,orders{ii});
        
        if ii==1
            
            ABCD=[];
            
            matsizes=[];
            
        end
        
    end
    
end

end
