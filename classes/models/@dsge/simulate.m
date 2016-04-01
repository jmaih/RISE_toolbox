function varargout=simulate(obj,varargin)

if isempty(obj)
    
    out=simulate@rise_generic(obj);
    
    [varargout{1:nargout}]=utils.miscellaneous.mergestructures(out,...
        struct('simul_sig',1,'simul_pruned',false,'simul_order',[]),...
        'simul_fbs_horizon',0);
    
    return
    
end

% simul_fbs_horizon: horizon of testing for violations when simulating

% simul_order is solution order to use for simulation. Will most likely
% create problems with optimal policy if the user tries to simulate an optimal policy model
% with order of approximation higher than 1  

[varargout{1:nargout}]=simulate@rise_generic(obj,varargin{:});

end