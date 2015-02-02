function Initcond=set_simulation_initial_conditions(obj)
% set_simulation_initial_conditions - sets the initial conditions for forecasting, simulation and irfs
%
% Syntax
% -------
% ::
%
%   Initcond=set_simulation_initial_conditions(obj)
%
% Inputs
% -------
%
% - **obj** [rise|dsge|svar|rfvar]: model object
%
% Outputs
% --------
%
% - **Initcond** [struct]: Initial conditions for simulation
%
% More About
% ------------
%
% Examples
% ---------
%
% See also: 



% complementarity: do it here so as to keep the memory light. Here we pick
% the second argument only!!!
%-----------------------------------------------------------
[~,complementarity]=complementarity_memoizer(obj);


simul_pruned=false;
simul_sig=0;
simul_order=1;
nlags=1;
k_future=0;
if isa(obj,'svar')
    nlags=obj.nlags;
elseif isa(obj,'dsge')
    k_future=max(obj.exogenous.shock_horizon);
    simul_sig=obj.options.simul_sig;
    simul_pruned=obj.options.simul_pruned;
    if isempty(obj.options.simul_order);
        simul_order=obj.options.solve_order;
    else
        simul_order=obj.options.simul_order;
    end
else
    error(['model of class ',class(obj),' not ready for simulation'])
end

% one initial condition despite multiple regimes
%------------------------------------------------
ss=cell2mat(obj.solution.ss);
[PAI,retcode]=initial_markov_distribution(obj.solution.transition_matrices.Q);
if retcode
    error(decipher(retcode))
end
ss=sum(bsxfun(@times,ss,PAI(:).'),2);
y0=struct('y',[],'y_lin',[]);
y0(1).y=vec(ss(:,ones(1,nlags)));

simul_history_end_date=0;
simul_historical_data=obj.options.simul_historical_data;
shocks=[];
states=[];
scale_shocks=1;
if ~isempty(simul_historical_data)
    if isstruct(simul_historical_data)
        simul_historical_data=ts.collect(simul_historical_data);
    end
    if ~isa(simul_historical_data,'ts')
        error('historical database must be a ts')
    end
    simul_history_end_date=obj.options.simul_history_end_date;
    if isempty(simul_history_end_date)
        simul_history_end_date=simul_historical_data.finish;
    end
    
    % set the endogenous variables
    %-----------------------------
     set_endogenous_variables()
        
    % now load shocks
    %----------------
    set_shocks_and_states()
    % in case shocks are not produced right above, set to 0 all the shocks
    % that will be created randomly below
    scale_shocks=0; 
end

Qfunc=prepare_transition_routine(obj);

if ~simul_pruned
    y0=rmfield(y0,'y_lin');
end
if ~isfield(obj.options,'simul_update_shocks_handle')
    obj.options.simul_update_shocks_handle=[];
end
if ~isfield(obj.options,'simul_do_update_shocks')
    obj.options.simul_do_update_shocks=[];
end

Initcond=struct('y',{y0},...
    'PAI',PAI,...
    'simul_history_end_date',simul_history_end_date,...
    'simul_sig',simul_sig,...
    'simul_order',simul_order,...
    'Qfunc',Qfunc,...
    'complementarity',complementarity,...
    'random',true,...
    'nsteps',obj.options.simul_periods,...
    'k_future',k_future,...
    'simul_update_shocks_handle',obj.options.simul_update_shocks_handle,...
    'simul_do_update_shocks',obj.options.simul_do_update_shocks,...
    'forecast_conditional_hypothesis',obj.options.forecast_conditional_hypothesis);
%-----------------------------------------
Initcond.burn=obj.options.simul_burn;
if ~isempty(shocks)
    % shocks have already been set from the initial conditions no
    % burn-in simulations necessary
    Initcond.burn=0;
else
    exo_nbr=sum(obj.exogenous.number);
    which_shocks=true(1,exo_nbr);
    which_shocks(obj.exogenous.is_observed)=false;
    shocks=scale_shocks*utils.forecast.create_shocks(exo_nbr,[],~which_shocks,Initcond);
end
if isempty(states)
    states=nan(Initcond.nsteps+Initcond.burn,1);
end
simul_regime=obj.options.simul_regime;
if all(isnan(states)) && ~isempty(simul_regime)
    nregs=numel(simul_regime);
    if nregs==1
        states(:)=simul_regime;
    else
        nregs=min(nregs,numel(states));
        states(1:nregs)=simul_regime(1:nregs);
        % extend the last regime
        %------------------------
        states(nregs+1:end)=states(nregs);
    end
end

%-----------------------------------------
Initcond.states=states;
Initcond.shocks=shocks;

    function set_shocks_and_states()
        % check there is no shock with name regime (this should be done right
        % from the parser)
        %--------------------------------------------------------------------
        shock_names=obj.exogenous.name;
        
        % for the shocks, the first state is right after the end of history
        %------------------------------------------------------------------
        shocks_start=utils.forecast.load_start_values(shock_names,simul_historical_data,...
            date2serial(simul_history_end_date)+1);
        regime_start=utils.forecast.load_start_values('regime',simul_historical_data,...
            date2serial(simul_history_end_date)+1);
        
        % put into a single page
        %-----------------------
        shocks_start=shocks_start(:,:);
        regime_start=regime_start(:,:);
        % extend the shocks with zeros
        %-----------------------------
        simul_periods=obj.options.simul_periods;
        nshocks=numel(shock_names)+1;
        shocks=zeros(nshocks,simul_periods+k_future);
        ncols=min(size(shocks_start,2),simul_periods+k_future);
        shocks(1:end-1,1:ncols)=shocks_start(:,1:ncols);
        ncols=min(size(regime_start,2),simul_periods+k_future);
        shocks(end,1:ncols)=regime_start(1,1:ncols);
        clear regime_start shocks_start
        states=shocks(end,:);
        shocks=shocks(1:end-1,:);
        % there is no such a thing as a zero state/regime
        states(states==0)=nan;
    end

    function set_endogenous_variables()
        endo_names=obj.endogenous.name(:);
        endo_names=endo_names(:,ones(1,nlags));
        for ilag=2:nlags
            endo_names(:,ilag)=strcat(endo_names(:,ilag),sprintf('{-%0.0f}',ilag-1));
        end
        endo_names=endo_names(:).';

        y0(1).y=utils.forecast.load_start_values(endo_names,simul_historical_data,...
            simul_history_end_date,y0(1).y);
    end
end