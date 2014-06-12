function [ovSolution,Q,PAI,State,Initcond,simul_burn,simul_periods,retcode]=...
    simulation_preamble(obj,varargin)
obj=set(obj,varargin{:});
if isa(obj,'dsge') && ...
        ~isempty(obj.options.simul_order) && ...
        obj.options.simul_order>obj.options.solve_order
    obj.options.solve_order=obj.options.simul_order;
end

% solution and system
%--------------------
[obj,retcode]=solve(obj);
if retcode
    if nargout==8
        ovSolution=[];
        Q=[];
        PAI=[];
        State=[];
        Initcond=[];
        simul_burn=[];
        simul_periods=[];
    else
        error([mfilename,':: no solution: the model cannot be simulated'])
    end
    return
end
% initial conditions
%-------------------
Initcond=generic_tools.set_simulation_initial_conditions(obj);

% load the solution matrices (Abstract method)
%--------------------------------------------
%[T,steady_state,state_vars_locations,y0]
ovSolution=load_order_var_solution(obj,Initcond.y);
simul_periods=obj.options.simul_periods;
simul_regime=obj.options.simul_regime;
if ~isempty(Initcond.shocks)
    % There is no burn-in period in this case
    %----------------------------------------
    obj.options.simul_burn=0;
    % also set the simulation of shocks to false
    %-------------------------------------------
    obj.options.simul_shocks=false;
    simul_regime=Initcond.shocks(end,:);
    % remove the regime row
    Initcond.shocks(end,:)=[];
end
simul_burn=obj.options.simul_burn;
Q=Initcond.Q;
PAI=Initcond.PAI;
% the following options are collected from the Initcond as they
% come from the dsge
%------------------------------------------------------------------------
if obj.options.simul_pruned
    ovSolution.y0.y_lin=[];
end

h=obj.markov_chains.regimes_number;
if ~isempty(simul_regime) && ~all(ismember(simul_regime,1:h))
    error([mfilename,':: all elements in the history of the markov states must be between 1 and ',int2str(h)])
end
State=generic_tools.set_simulation_regimes(simul_regime,simul_periods,simul_burn);

Initcond=rmfield(Initcond,{'PAI','Q','y'});

end
