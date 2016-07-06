function obj=estimate(obj,varargin)

if isempty(obj)
    
    obj=estimate@generic(obj);
    
    obj.estim_scale_free=true;
    
    obj.estim_no_priors=false;
    
    return
    
end

obj=set(obj,varargin{:});

obj=data_location_dispatch(obj);

obj=vstar.set_baseline_parameters(obj);

obj=estimate@generic(obj);

end
