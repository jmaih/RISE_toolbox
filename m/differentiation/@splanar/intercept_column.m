function obj=intercept_column(obj,pointer)
% intercept_column - builds a scalar splanar object from a
% vectorized splanar object. The pointer argument points the
% element in the vector to be used.

if isnumeric(obj) && numel(obj.func)>1
    
    obj.func=obj.func(pointer);
    
elseif ~isempty(obj.args)
    
    is_simplifiable=false;
    
    incid=sparse(1,numel(obj.incidence));
    
    for iarg=1:numel(obj.args)
        
        if ~isa(obj.args{iarg},'splanar')
            
            continue
            
        end
        
        obj.args{iarg}=intercept_column(obj.args{iarg},pointer);
        
        if ~is_simplifiable
            
            argifunc=obj.args{iarg}.func;
            
            is_simplifiable=isnumeric(argifunc)&&...
                isscalar(argifunc) && any(abs(argifunc)-[0,1]==0);
            
        end
        
        if ~isempty(obj.args{iarg}.incidence)
            
            incid=incid|obj.args{iarg}.incidence;
            
        end
        
    end
    % the line below will correct for zeros, ones, etc. as well as
    % incidences.
    if is_simplifiable
        
        obj=feval(obj.func,obj.args{:});
        
    else
        
        obj.incidence=incid;
        
    end
    
end

end
