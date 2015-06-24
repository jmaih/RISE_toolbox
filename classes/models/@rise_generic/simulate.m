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
%
%   - **simul_periods** [integer|{100}]: number of simulation periods
%
%   - **simul_burn** [integer|{100}]: number of burn-in periods
%
%   - **simul_algo** [[{mt19937ar}| mcg16807|mlfg6331_64|mrg32k3a|
%       shr3cong|swb2712]]: matlab's seeding algorithms
%
%   - **simul_seed** [numeric|{0}]: seed of the computations
%
%   - **simul_historical_data** [ts|struct|{''}]: historical data from
%       which the simulations are based. If empty, the simulations start at
%       the steady state.
%
%   - **simul_history_end_date** [char|integer|serial date]: last date of
%       history
%
%   - **simul_regime** [integer|vector|{[]}]: regimes for which the model
%       is simulated
%
%   - **simul_update_shocks_handle** [function handle]: we may want to
%       update the shocks if some condition on the state of the economy is
%       satisfied. For instance, shock monetary policy to keep the interest
%       rate at the floor for an extented period of time if we already are
%       at the ZLB/ZIF. simul_update_shocks_handle takes as inputs the
%       current shocks and the state vector (all the endogenous variables)
%       and returns the updated shocks. But for all this to be put into
%       motion, the user also has to turn on **simul_do_update_shocks** by
%       setting it to true.
%
%   - **simul_do_update_shocks** [true|{false}]: update the shocks based on
%       **simul_update_shocks_handle** or not. 
%
%   - **simul_to_time_series** [{true}|false]: if true, the output is a
%       time series, else a cell array with a matrix and information on
%       elements that help reconstruct the time series.
%
%   - **simul_honor_constraints** [true|{false}]: honor restrictions during
%       simulations. If true, agents have to be able anticipate the future.
%
%   - **simul_frwrd_back_shoot** [true|{false}]: uses an algorithm that
%   checks the constraints are satisfied at all horizons instead of just
%   one at a time.
%
%   - **simul_shock_uncertainty** [{true}|false]: draw shocks over the
%   simulation horizon.
%
%   - **simul_honor_constraints_through_switch** [true|{false}]: if true,
%   constraints are honored through the switching mechanism. In that case
%   the number of regimes should be greater than 1. If false, constraints
%   are honored through an anticipatory behavior. In this case, there
%   should be shocks that are foreseen.
%
%   - **simul_anticipate_zero** [true|{false}]: When shocks are drawn, this
%   option allows to impose that agents continue to see only the
%   contemporaneous shocks.
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
        'simul_honor_constraints',false,...
        'simul_frwrd_back_shoot',false,...
        'simul_to_time_series',true,...
        'simul_shock_uncertainty',true,...
        'simul_honor_constraints_through_switch',false,...
        'simul_anticipate_zero',false);
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
do_dsge_var=false;
if isa(obj,'dsge')
    do_dsge_var=obj.is_dsge_var_model && obj.options.dsgevar_var_regime;
end
[T,~,steady_state,new_order,state_vars_location]=load_solution(obj,'ov',do_dsge_var);

Initcond=utils.forecast.initial_conditions_to_order_var(Initcond,new_order,obj.options);

y0=Initcond.y;
if numel(y0)>1
    error('more than one initial conditions')
end

[y,states,retcode,~,myshocks]=utils.forecast.multi_step(y0(1),steady_state,T,...
    state_vars_location,Initcond);

% add initial conditions: (only the actual data on the first page)
% ONLY IF THERE WAS NO BURN-IN
%---------------------------------------------------------------
if Initcond.burn==0
    % we include history
    y=[y0(1).y(:,:,1),y];
    y0cols=size(y0(1).y,2);
else
    % then we start after the end of history
    y0cols=0;
end

% recompose the output if we have a dsge-var
%-------------------------------------------
if do_dsge_var
    endo_nbr=obj.endogenous.number;
    max_rows=obj.observables.number(1);
    relevant=obj.inv_order_var(obj.observables.state_id);
    ncols=size(y,2);
    yold=y;
    y=zeros(endo_nbr,ncols)+1i;
    y(relevant,:)=yold(1:max_rows,:);
end

% put y in the correct order before storing
%------------------------------------------
y=re_order_output_rows(obj,y);

start_date=serial2date(date2serial(Initcond.simul_history_end_date)-y0cols+1);

smpl_states=numel(states);
[states_,markov_chains]=regimes2states(states);
smply=size(y,2);
[exo_nbr,smplx]=size(myshocks);
myshocks=[zeros(exo_nbr,y0cols),myshocks];
smplx=smplx+y0cols;
vnames=[obj.endogenous.name,obj.exogenous.name,'regime',markov_chains];
endo_nbr=obj.endogenous.number;
smpl=max(smplx,smply);
yy=nan(smpl,endo_nbr+exo_nbr+1+numel(markov_chains));
yy(1:smply,1:endo_nbr)=y(1:endo_nbr,:)';
yy(1:smplx,endo_nbr+(1:exo_nbr))=myshocks';
yy(y0cols+1:smply,endo_nbr+exo_nbr+1:end)=[states,states_];
if obj.options.simul_to_time_series
    % store the simulations in a database: use the date for last observation in
    % history and not first date of forecast
    %--------------------------------------------------------------------------
    % store only the relevant rows in case we are dealing with a VAR with many
    % lags
    db=ts(start_date,yy,vnames);
    db=pages2struct(db);
else
    db={yy,...
        {
        '1=time, start_date=',start_date
        '2= endogenous names',vnames
        }
        };
end

    function [states,markov_chains]=regimes2states(regimes_history)
        regimes_tables=obj.markov_chains.regimes;
        markov_chains=regimes_tables(1,2:end);
        nchains=numel(markov_chains);
        reg_table=cell2mat(regimes_tables(2:end,2:end));
        states=nan(smpl_states,nchains);
        
        max_reg=max(regimes_history);
        for ireg=1:max_reg
            pos=regimes_history==ireg;
            n=sum(pos);
            states(pos,:)=reg_table(ireg*ones(n,1),:);
        end
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