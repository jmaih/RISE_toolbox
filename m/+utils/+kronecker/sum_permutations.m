function P=sum_permutations(ABCD,matsizes,options,varargin)
% sum_permutations -- sums the permutations of a tensor product
%
% ::
%
%
%   P=sum_permutations(ABCD,matsizes,[],order_1,order_2,...,order_n)
%
%   P=sum_permutations(ABCD,matsizes,[],{order_1,order_2,...,order_n})
%
%   P=sum_permutations(ABCD,matsizes,options,order_1,order_2,...,order_n)
%
%   P=sum_permutations(ABCD,matsizes,options,{order_1,order_2,...,order_n})
%
% Args:
%
%    - **ABCD** [matrix]: tensor product of k matrices
%
%    - **matsizes** [matrix]: k x 2 matrix of size of the various matrices
%    entering the tensor. Each row represents the size of a matrix and it is
%    assumed that the main(or first) tensor product is ordered [1,2,...,k]
%
%    - **options** [empty|struct]: structure with various options such as
%      - **algo** [{'shuf1'}|'shuf2'|'old'|'perm']: shuf1 shuffles the
%      matrix, shuf2 pre and post multiplies the matrix, old may or may not
%      construct a grid, perm uses a permutation-type of strategy.
%      - **use_grid** [true|{false}]: use grid in the old algorithm:a grid is
%          used to compute the indexes of the main kronecker product
%      - **skip_first** [true|{false}]: if true, the original input matrix is
%      not added to the sum
%
%    - **order_i** [vector]: permutation of [1,2,...,k]. N.B: all elements
%    1,2,...,k should be part of the vector and the vector should have exactly
%    k elements
%
% Returns:
%    :
%
%    - **P** sum of ABCD and its permutations
%
% Note:
%
% Example:
%
%    See also: tensorperm

default_options={
    'algo','shuf1',@(x)ismember(x,{'shuf1','shuf2','old','perm'}),'algo must be "shuf1", "shuf2", "old" or "perm"'
    'use_grid',false,@(x)islogical(x),'use_grid must be a logical'
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

switch options.algo
    
    case 'old'
        
        old_algorithm()
        
    case 'perm'
        
        permute_algorithm()
        
    case 'shuf1'
        
        shuffle1_algorithm()
        
    case 'shuf2'
        
        shuffle2_algorithm()
        
end

    function shuffle1_algorithm()
        
        for ii=1:numel(orders)
            
            P=P+utils.kronecker.shuffle_tensor1(ABCD,matsizes,orders{ii});
            
        end
        
    end

    function shuffle2_algorithm()
        
        for ii=1:numel(orders)
            
            P=P+utils.kronecker.shuffle_tensor2(ABCD,matsizes,orders{ii});
            
        end
        
    end

    function permute_algorithm()
        
        for ii=1:numel(orders)
            
            P=P+utils.kronecker.permute_tensor(ABCD,matsizes,orders{ii});
            
            if ii==1
                
                ABCD=[];
                
                matsizes=[];
                
            end
            
        end
        
    end

    function old_algorithm()
        
        if options.use_grid
            
            [irows,jcols]=utils.kronecker.tensorperm(matsizes,orders{:},'grid');
            
        else
            
            [irows,jcols]=utils.kronecker.tensorperm(matsizes,orders{:});
            
        end
        
        for ii=1:size(irows,2)
            
            P=P+ABCD(irows(:,ii),jcols(:,ii));
            
        end
        
    end

end
