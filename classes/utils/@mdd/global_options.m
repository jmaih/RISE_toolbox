function opts=global_options(opts0)

if nargin==0||isempty(opts0)
    
    opts0=struct();
    
end

defaults0={ % arg_names -- defaults -- checks -- error_msg
    'center_at_mean',false,@(x)isscalar(x) && islogical(x),...
    'center_at_mean should be a logical scalar'
    
    'L',500,@(x)isnumeric(x) && isreal(x) && isfinite(x) && ...
    ceil(x)==floor(x) && x>0,...
    'L (# of IID draws) should be an integer'
    
    'debug',false,@(x)isscalar(x) && islogical(x),'debug should be a logical scalar'
    
    'draws_are_iid',false,@(x)isscalar(x) && islogical(x),...
    'draws_are_iid should be a logical scalar'
    
    };

defaults1=fix_point_iterator();

defaults1(:,1)=strrep(defaults1(:,1),'(r)','');

defaults=[defaults0;defaults1];

opts=cell2struct(defaults(:,2),defaults(:,1),1);

fields=fieldnames(opts0);

for ifield=1:numel(fields)
    
    f=fields{ifield};
    
    if ~isfield(opts,f)
        
        error(['"',f,'" is not a valid option for class mdd'])
        
    end
    
    fval=opts0.(f);
    
    if isempty(fval)
        
        continue
        
    end
    
    pos=strcmp(f,defaults(:,1));
    
    test=defaults{pos,3};
    
    answer=defaults{pos,4};
    
    assert(test(fval),answer)
    
    opts.(f)=fval;
    
end


end
