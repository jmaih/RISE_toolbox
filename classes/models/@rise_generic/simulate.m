function [db,State] = simulate(obj,varargin)
% simul_shocks can be:
%              - true : simul_periods and simul_burn are active
%              - false : simul_periods is active, simul_burn=0
% simul_historical_data contains the historical data as well as information
%              over the forecast horizon. This also includes the future
%              path of the regimes, which is denoted by "regime" in the
%              database.
% simul_history_end_date controls the end of history

% the random numbers specifications and algorithms are specified outside
% this function. Need to check this again!!!
if isempty(obj)
    if nargout>1
        error([mfilename,':: when the object is emtpy, nargout must be at most 1'])
    end
    db=struct('simul_periods',100,...
        'simul_burn',100,...
        'simul_algo','mt19937ar',...  % [{mt19937ar}| mcg16807| mlfg6331_64|mrg32k3a|shr3cong|swb2712]
        'simul_seed',0,...
        'simul_historical_data',ts.empty(0),...
        'simul_history_end_date','',...
        'simul_shocks',true,... [true,false]
        'simul_start_date','',...
        'simul_regime',[]);
    return
end
nobj=numel(obj);
if nobj>1
    State=cell(1,nobj);
    db=cell(1,nobj);
    for iobj=1:nobj
        [db{iobj},State{iobj}] = simulate(obj(iobj),varargin{:});
    end
    db=format_simulated_data_output(db);
    return
end

        [obj,retcode]=solve(obj);
        % note that the solving of the model may change the perturbation
        % order. More explicitly optimal policy irfs will be computed for a
        % perturbation of order 1 no matter what order the user chooses.
        % This is because we have not solved the optimal policy problem
        % beyond the first order.
        if retcode
            error('model cannot be solved')
        end

        % initial conditions
        %-------------------
        Initcond=generic_tools.set_simulation_initial_conditions(obj);
        PAI=Initcond.PAI;
        Q=Initcond.Q;

        ovSolution=load_order_var_solution(obj,Initcond.y);
        y0=ovSolution.y0;
        T=ovSolution.T;
        steady_state=ovSolution.steady_state;
        
        % here we need to start at one single point: and so we aggregate y0
        %------------------------------------------------------------------
        [y0]=utils.forecast.aggregate_initial_conditions(PAI,y0);
        
        % further options
        %----------------
        simul_order=obj.options.simul_order;
        if isempty(simul_order)
            simul_order=obj.options.solve_order;
        end
        options=struct('simul_sig',obj.options.simul_sig,...
            'simul_order',simul_order,...
            'burn',obj.options.simul_burn,...
            'nsimul',1,...
            'impulse',[],...
            'random',true,...
            'girf',false,...
            'nsteps',obj.options.simul_periods);
        if isa(obj,'dsge')
            options.k_future=max(obj.exogenous.shock_horizon);
        else
            options.k_future=0;
        end
        
        shocks=Initcond.shocks;
        if ~isempty(shocks)
            options.burn=0;
        else
            exo_nbr=sum(obj.exogenous.number);
            which_shocks=true(1,exo_nbr);
            which_shocks(obj.exogenous.is_observed)=false;
            shocks=utils.forecast.create_shocks(exo_nbr,[],~which_shocks,options);
        end
        states=nan(obj.options.simul_periods,1);
        [y,states,retcode]=utils.forecast.multi_step(y0,steady_state,T,shocks,states,Q,PAI,options);

        % initialize output
        %------------------

%       State=State(simul_burn+1:end);

% put y in the correct order before storing
%------------------------------------------
if isa(obj,'dsge')
    y=y(obj.inv_order_var.after_solve,:);
end

% store the simulations in a database: use the date for last observation in
% history and not first date of forecast
%--------------------------------------------------------------------------
y0cols=size(y0(1).y,2);
start_date=serial2date(date2serial(0)-y0cols+1);
db=ts(start_date,y',obj.endogenous.name);
db=pages2struct(db);

end

function sim_data=format_simulated_data_output(sim_data)
nobj=numel(sim_data);
if nobj==1
    sim_data=sim_data{1};
else
    if isempty(sim_data{1})
        return
    end
    sim_data=utils.time_series.concatenate_series_from_different_models(sim_data);
end
end

% function [db,State] = simulate(obj,varargin)
% % simul_shocks can be:
% %              - true : simul_periods and simul_burn are active
% %              - false : simul_periods is active, simul_burn=0
% % simul_historical_data contains the historical data as well as information
% %              over the forecast horizon. This also includes the future
% %              path of the regimes, which is denoted by "regime" in the
% %              database.
% % simul_history_end_date controls the end of history
% 
% % the random numbers specifications and algorithms are specified outside
% % this function. Need to check this again!!!
% if isempty(obj)
%     if nargout>1
%         error([mfilename,':: when the object is emtpy, nargout must be at most 1'])
%     end
%     db=struct('simul_periods',100,...
%         'simul_burn',100,...
%         'simul_algo','mt19937ar',...  % [{mt19937ar}| mcg16807| mlfg6331_64|mrg32k3a|shr3cong|swb2712]
%         'simul_seed',0,...
%         'simul_historical_data',ts.empty(0),...
%         'simul_history_end_date','',...
%         'simul_shocks',true,... [true,false]
%         'simul_start_date','',...
%         'simul_regime',[]);
%     return
% end
% nobj=numel(obj);
% if nobj>1
%     State=cell(1,nobj);
%     db=cell(1,nobj);
%     for iobj=1:nobj
%         [db{iobj},State{iobj}] = simulate(obj(iobj),varargin{:});
%     end
%     db=format_simulated_data_output(db);
%     return
% end
% 
% [ovSolution,Q,PAI,State,Initcond,simul_burn,simul_periods]=simulation_preamble(obj,varargin{:});
% 
% y0=ovSolution.y0;
% T=ovSolution.T;
% steady_state=ovSolution.steady_state;
% state_vars_location=ovSolution.state_vars_location;
% clear ovSolution
% 
% y0cols=size(y0.y,2);
% random=isequal(obj.options.simul_shocks,true) && isempty(obj.options.simul_historical_data);
% nx=sum(obj.exogenous.number);
% shock_id=[];
% impulse=[];
% k=max(obj.exogenous.shock_horizon);
% det_vars=obj.exogenous.is_observed;
% absent_shocks=isempty(Initcond.shocks);
% endo_nbr=obj.endogenous.number(end);
% for t=1:simul_burn+simul_periods
%     % draw a state
%     %-------------
%     [State(t+1),Q,PAI,retcode]=generic_tools.choose_state(State(t+1),Q,...
%         PAI,y0.y(1:endo_nbr,end));
%     if retcode
%         if obj.options.debug
%             utils.error.decipher(retcode)
%         end
%         return
%     end
%     if absent_shocks
%         shocks_t=utils.forecast.create_shocks(nx,random,shock_id,impulse,k,det_vars);
%     else
%         shocks_t=Initcond.shocks(:,t+(0:k));
%     end
%     % Simulate one step
%     %------------------
%     rt=State(t+1);
%     y1=utils.forecast.one_step(T(:,rt),...
%         y0,...
%         steady_state{rt},...
%         state_vars_location,...
%         Initcond.simul_sig,...
%         shocks_t,...
%         Initcond.simul_order);
%     
%     if t>simul_burn
%         if t==simul_burn+1
%             y=[y0.y(1:endo_nbr,:),zeros(endo_nbr,simul_periods)];
%         end
%         y(:,t-simul_burn+y0cols)=y1.y(1:endo_nbr,end);
%     end
%     y0=y1;
% end
% State=State(simul_burn+1:end);
% 
% % put y in the correct order before storing
% %------------------------------------------
% if isa(obj,'dsge')
%     y=y(obj.inv_order_var.after_solve,:);
% end
% 
% % store the simulations in a database: use the date for last observation in
% % history and not first date of forecast
% %--------------------------------------------------------------------------
% start_date=serial2date(date2serial(Initcond.simul_history_end_date)-y0cols+1);
% db=ts(start_date,y',obj.endogenous.name);
% db=pages2struct(db);
% 
% end
% 
% 
% function sim_data=format_simulated_data_output(sim_data)
% nobj=numel(sim_data);
% if nobj==1
%     sim_data=sim_data{1};
% else
%     if isempty(sim_data{1})
%         return
%     end
%     sim_data=utils.time_series.concatenate_series_from_different_models(sim_data);
% end
% end