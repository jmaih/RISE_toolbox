function out=A_times_reordered_kron_Q1_Qk(A,matsizes,orders,options,varargin)
% matsizes : sizes of matrices entering kron(Q1,...,Qk)
% orders : orders of permutations of the Q1,...,Qk matrices

default_options={
    'skip_first',false,@(x)islogical(x),'skip_first must be a logical'
    };

if isempty(options)
    
    options=cell2struct(default_options(:,2),default_options(:,1),1);
    
else
    
    if ~isstruct(options)
        
        error('options must be a structure or empty')
        
    end
    
    options=parse_arguments(default_options,options);
    
end

engine=@utils.kronecker.A_times_kron_Q1_Qk_master;

nmat=length(varargin);

siza=zeros(nmat,2);

for ii=1:nmat
    
    siza(ii,:)=size(varargin{ii});
    
end

siza=prod(siza,1);

if options.skip_first
    
    out=0;
    
else
    
    out=engine('fast',A,varargin{:});
    
end

for ii=1:numel(orders)
    
    newOrder=orders{ii};
    
    [lm,rm]=utils.kronecker.find_reordering(siza,matsizes,newOrder);
    
    tmp=engine('fast',A*lm,varargin{:});
    
    out=out+tmp*rm;
    
end

end