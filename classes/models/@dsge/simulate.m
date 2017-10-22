function varargout=simulate(obj,varargin)

if isempty(obj)
    
    mydefaults=the_defaults();
    
    mydefaults=[mydefaults
        simulate@generic(obj)];
    
    if nargout
        
        varargout={mydefaults};
        
    else
        
        disp_defaults(mydefaults);
        
    end
    
    return
    
end

% simul_fbs_horizon: horizon of testing for violations when simulating

% simul_order is solution order to use for simulation. Will most likely
% create problems with optimal policy if the user tries to simulate an optimal policy model
% with order of approximation higher than 1  

[varargout{1:nargout}]=simulate@generic(obj,varargin{:});

end

function d=the_defaults()

num_fin=@(x)isnumeric(x) && isscalar(x) && isfinite(x);

num_fin_int=@(x)num_fin(x) && floor(x)==ceil(x) && x>=0;

d={
    'simul_sig',1,@(x)num_fin_int(x) && x>=0,'simul_sig must be >=0'
    
    'simul_pruned',false,@(x)islogical(x),'simul_pruned must be a logical'
    
    'simul_order',[],@(x)num_fin_int(x) && x>=1,'simul_order must be a positive integer'
    
    'simul_fbs_horizon',0,@(x)num_fin_int(x),'simul_fbs_horizon must be a finite and positive integer'
    };

end