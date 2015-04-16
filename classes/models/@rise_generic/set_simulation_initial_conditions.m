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
% - if future values of endogenous are found, they are discarded the
% variables are not explicitly declared as conditioning endogenous
% variables.
%
% - if either endogenous or exogenous variables are declared as
% conditioning variables, we automatically have a conditional forecasting
% exercise.
%
% - if no variable is declared as conditioning, we have a simulation
% exercise, in which case, there is a burn-in period.
%
% - Initialization is always done at the steady state and the declared
% endogenous variables merely override the steady state values. This means
% that if a variable is not found in the database, it is initialized at its
% steady state.
%
% Examples
% ---------
%
% See also: 



% complementarity: do it here so as to keep the memory light
%-----------------------------------------------------------
[sep_compl,complementarity]=complementarity_memoizer(obj);

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

exo_nbr=sum(obj.exogenous.number);
simul_history_end_date=0;
simul_historical_data=obj.options.simul_historical_data;
shocks=[];
states=[];

is_conditional_forecasting=false;
if ~isempty(simul_historical_data)
    % no burn-in with historical data
    %----------------------------------
    obj.options.simul_burn=0;
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

% shock structure: initial+anticipate
shock_structure=false(exo_nbr,k_future);
if k_future
    for iexo=1:exo_nbr
        shock_structure(iexo,1:obj.exogenous.shock_horizon(iexo))=true;
    end
end
shock_structure=[false(exo_nbr,1),shock_structure];
Initcond=struct('y',{y0},...
    'PAI',PAI,...
    'simul_history_end_date',simul_history_end_date,...
    'simul_sig',simul_sig,...
    'simul_order',simul_order,...
    'Qfunc',Qfunc,...
    'complementarity',complementarity,...
    'random',obj.options.simul_shock_uncertainty,...
    'nsteps',obj.options.simul_periods,...
    'k_future',k_future,...
    'simul_update_shocks_handle',obj.options.simul_update_shocks_handle,...
    'simul_do_update_shocks',obj.options.simul_do_update_shocks,...
    'forecast_conditional_hypothesis',obj.options.forecast_conditional_hypothesis,...
    'shock_structure',shock_structure,...
    'sep_compl',sep_compl);
%-----------------------------------------
Initcond.burn=obj.options.simul_burn;
if is_conditional_forecasting
    % shocks have already been set from the initial conditions no
    % burn-in simulations necessary
    Initcond.burn=0;
else
    which_shocks=true(1,exo_nbr);
    which_shocks(obj.exogenous.is_observed)=false;
    shocks=utils.forecast.create_shocks(exo_nbr,[],~which_shocks,Initcond);
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
        shock_names=obj.exogenous.name;
        conditional_shocks=obj.options.forecast_cond_exo_vars;
        simul_periods=obj.options.simul_periods;
        states=nan(1,simul_periods);
        shocks=nan(exo_nbr,simul_periods+k_future);
        if isempty(conditional_shocks)
            is_conditional_forecasting=size(y0(1).y,3)>1;
        else
            is_conditional_forecasting=true;
            if ischar(conditional_shocks)
                conditional_shocks=cellstr(conditional_shocks);
            end
            regime_pos=find(strcmp('regime',conditional_shocks));
            if ~isempty(regime_pos)
                conditional_shocks(regime_pos)=[];
            end
            shock_pos=[];
            if ~isempty(conditional_shocks)
                shock_pos=locate_variables(conditional_shocks,shock_names);
            end
            if ~isempty(regime_pos)
                conditional_shocks=[conditional_shocks(:).','regime'];
            end
            % the length/number of pages of the dataset depends on the
            % horizon of the shocks but for the moment, we consider all pages
            pages=[];
            verdier=utils.time_series.data_request(simul_historical_data,...
                conditional_shocks,date2serial(simul_history_end_date),...
                [],pages);
            % remove first page (shocks of the past!!!)
            %-------------------------------------------
            % squeeze does not work well with singleton dimensions...
            verdier=permute(verdier(:,:,2:end),[1,3,2]); % <-- squeeze(verdier(:,:,2:end));
            npages=size(verdier,2);
            % now chop until the forecast horizon
            %--------------------------------------
            ncols=min(simul_periods,npages);
            if ~isempty(regime_pos)
                states(1:ncols)=verdier(end,1:ncols);
                verdier(end,:)=[];
            end
            ncols=min(simul_periods+k_future,npages);
            shocks(shock_pos,1:ncols)=verdier(:,1:ncols);
        end
        if ~is_conditional_forecasting
            states=[];
            shocks=[];
        end
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
        
        % past regimes for the computation of initial regime probabilities
        %-----------------------------------------------------------------
        PAI_lag=utils.forecast.load_start_values(obj.markov_chains.regime_names,...
            simul_historical_data,simul_history_end_date,nan(size(PAI)));
        if all(~isnan(PAI_lag))
            PAI=transpose(obj.solution.transition_matrices.Q)*PAI_lag;
        end

        conditional_vars=obj.options.forecast_cond_endo_vars;
        % x it future values of variables not declared as conditioning
        %------------------------------------------------------------------
        if size(y0(1).y,3)>1
            is_cond_var=ismember(obj.endogenous.name,conditional_vars);
            y0(1).y(~is_cond_var,:,2:end)=nan;
            % start trimming from the back
            %-----------------------------
            while size(y0(1).y,3)>1
                last=y0(1).y(:,:,end);
                if any(~isnan(last(:)))
                    break
                end
                y0(1).y(:,:,end)=[];
            end
        end
    end
end