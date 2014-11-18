function [db,states,retcode] = simulate(obj,varargin)
% simulate - simulates a RISE model
%
% Syntax
% -------
% ::
%
%   [db,states,retcode] = simulate(obj,varargin)
%
% Inputs
% -------
%
% - **obj** [rfvar|dsge|rise|svar]: model object
%
% - **varargin** : additional arguments including but not restricted to
%   - **simul_periods** [integer|{100}]: number of simulation periods
%   - **simul_burn** [integer|{100}]: number of burn-in periods
%   - **simul_algo** [[{mt19937ar}| mcg16807|mlfg6331_64|mrg32k3a|
%       shr3cong|swb2712]]: matlab's seeding algorithms
%   - **simul_seed** [numeric|{0}]: seed of the computations
%   - **simul_historical_data** [ts|struct|{''}]: historical data from
%       which the simulations are based. If empty, the simulations start at
%       the steady state.
%   - **simul_history_end_date** [char|integer|serial date]: last date of
%       history
%   - **simul_regime** [integer|vector|{[]}]: regimes for which the model
%       is simulated
%   - **simul_update_shocks_handle** [function handle]: we may want to
%       update the shocks if some condition on the state of the economy is
%       satisfied. For instance, shock monetary policy to keep the interest
%       rate at the floor for an extented period of time if we already are
%       at the ZLB/ZIF. simul_update_shocks_handle takes as inputs the
%       current shocks and the state vector (all the endogenous variables).
%       The user also has to turn on **simul_do_update_shocks** by setting
%       it to true.
%   - **simul_do_update_shocks** [true|{false}]: update the shocks based on
%       **simul_update_shocks_handle** or not. 
%   - **simul_to_time_series** [{true}|false]: if true, the output is a
%       time series, else a cell array with a matrix and information on
%       elements that help reconstruct the time series.
%
% Outputs
% --------
%
% - **db** [struct|cell array]: if **simul_to_time_series** is true, the
%   output is a time series, else a cell array with a matrix and
%   information on elements that help reconstruct the time series.
%
% - **states** [vector]: history of the regimes over the forecast horizon
%
% - **retcode** [integer]: if 0, the simulation went fine. Else something
%   got wrong. In that case one can understand the problem by running
%   decipher(retcode)
%
% More About
% ------------
%
% - **simul_historical_data** contains the historical data as well as
%   conditional information over the forecast horizon. It may also include
%   as an alternative to **simul_regime**, a time series with name
%   **regime**, which indicates the regimes over the forecast horizon.
%
% Examples
% ---------
%
% See also: 


if isempty(obj)
    if nargout>1
        error([mfilename,':: when the object is emtpy, nargout must be at most 1'])
    end
    db=struct('simul_periods',100,...
        'simul_burn',100,...
        'simul_algo','mt19937ar',... 
        'simul_seed',0,...
        'simul_historical_data',ts.empty(0),...
        'simul_history_end_date','',...
        'simul_regime',[],...
        'simul_update_shocks_handle',[],...
        'simul_do_update_shocks',false,...
        'simul_to_time_series',true);
%         'simul_start_date','',... does not seem to be in use
    return
end
nobj=numel(obj);
if nobj>1
    retcode=nan(1,nobj);
    states=cell(1,nobj);
    db=cell(1,nobj);
    for iobj=1:nobj
        [db{iobj},states{iobj},retcode(iobj)] = simulate(obj(iobj),varargin{:});
    end
    db=format_simulated_data_output(db);
    return
end

[obj,retcode]=solve(obj,varargin{:});
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
Initcond=set_simulation_initial_conditions(obj);

% load the order_var solution
%-----------------------------
[T,~,steady_state,new_order,state_vars_location]=load_solution(obj,'ov');
y0=Initcond.y;
if numel(y0)>1
    error('more than one initial conditions')
end
y0.y=full(y0.y); % [variables,ncols,1+n_conditions]
y0.y=y0.y(new_order,:,:); % [variables,ncols,1+n_conditions]
% adjust the transition and complementarity functions according to the
% order_var: they will be fed with order_var, but they expect the
% alphabetical order and so the order_var they are fed with has to be
% inverted prior to evaluation.
%--------------------------------------------------------------------------
iov(new_order)=1:numel(new_order);
Initcond.Qfunc=@(x)Initcond.Qfunc(x(iov));
Initcond.complementarity=@(x)Initcond.complementarity(x(iov));

[y,states,retcode]=utils.forecast.multi_step(y0(1),steady_state,T,...
    state_vars_location,Initcond);

% add initial conditions: only the actual data on the first page
%---------------------------------------------------------------
y=[y0(1).y(:,:,1),y];

% put y in the correct order before storing
%------------------------------------------
y=re_order_output_rows(obj,y);

y0cols=size(y0(1).y,2);
start_date=serial2date(date2serial(Initcond.simul_history_end_date)-y0cols+1);
if obj.options.simul_to_time_series
    % store the simulations in a database: use the date for last observation in
    % history and not first date of forecast
    %--------------------------------------------------------------------------
    % store only the relevant rows in case we are dealing with a VAR with many
    % lags
    db=ts(start_date,y(1:obj.endogenous.number(end),:)',obj.endogenous.name);
    db=pages2struct(db);
else
    db={y(1:obj.endogenous.number(end),:)',...
        {
        '1=time, start_date=',start_date
        '2= endogenous names',obj.endogenous.name
        }
        };
end
end

function sim_data=format_simulated_data_output(sim_data)
nobj=numel(sim_data);
if nobj==1
    sim_data=sim_data{1};
else
    if isempty(sim_data{1})||iscell(sim_data{1})
        return
    end
    sim_data=utils.time_series.concatenate_series_from_different_models(sim_data);
end
end