function final_names=create_variable_names(n,prefix,names0)

if nargin<3
    
    names0=[];
    
end

if numel(names0)>n
    
    error('more variables than required')
    
elseif numel(names0)==n
    
    final_names=names0;
    
else
    
    final_names=parser.create_state_list(prefix,n);
        
    if ~isempty(names0)
        
        final_names(1:numel(names0))=names0;
        
    end
    
end

end