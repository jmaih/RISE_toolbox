function F=fisher(obj,varargin) 

if isempty(obj)

    F=cell(0,4);

    return

end

obj.is_fisher=true;

if isempty(obj.linear_restrictions_data)
    
    obj=setup_restrictions(obj);
    
end

[~,F]=hessian(obj,varargin{:});

end