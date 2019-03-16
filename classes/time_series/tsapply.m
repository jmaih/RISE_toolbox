function db=tsapply(db,func)

% applies a basic operation to a ts object or a structure of ts objects.
% e.g.
% - db=tsapply(db,'log')
% - db=tsapply(db,@log)
% - db=tsapply(db,@(x)100*x)

if isstruct(db)
    
    fields=fieldnames(db);
    
    for ii=1:numel(fields)
        
        v=fields{ii};
        
        if isstruct(db.(v))
            
            db.(v)=tsapply(db.(v),func);
            
        elseif isa(db.(v),'ts')
            
            db.(v)=apply(db.(v),func);
            
        else
            
            error(['field ',v,' should be a ts object or a structure'])
            
        end
        
    end
    
elseif isa(db,'ts')
    
    db=apply(db,func);
    
else
    
    error('first input must be either a ts object or a structure of ts objects')
    
end

end