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
do_dsge_var=false;
if isa(obj,'svar')
    nlags=obj.nlags;
elseif isa(obj,'dsge')
    k_future=max(obj.exogenous.shock_horizon(:));
    simul_sig=obj.options.simul_sig;
    simul_pruned=obj.options.simul_pruned;
    if isempty(obj.options.simul_order);
        simul_order=obj.options.solve_order;
    else
        simul_order=obj.options.simul_order;
    end
    do_dsge_var=obj.is_dsge_var_model && obj.options.dsgevar_var_regime;
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
if do_dsge_var
    nlags=obj.options.dsgevar_lag;
    ss=ss(obj.observables.state_id,1);
end
y0=struct('y',[],'y_lin',[],'ycond',[],'econd',[]);
y0(1).y=vec(ss(:,ones(1,nlags)));

exo_nbr=sum(obj.exogenous.number);
simul_history_end_date=0;
simul_historical_data=obj.options.simul_historical_data;

is_conditional_forecasting=false;
has_data=~isempty(simul_historical_data);
if has_data
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
end

% set the endogenous variables
%-----------------------------
set_endogenous_variables()

% now load shocks
%----------------
set_shocks_and_states()

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
		maxlength=max(obj.exogenous.shock_horizon(:,iexo));
        shock_structure(iexo,1:maxlength)=true;
    end
end
shock_structure=[false(exo_nbr,1),shock_structure];

% load the options for utils.forecast.rscond.forecast
Initcond=utils.miscellaneous.reselect_options(obj.options,@utils.forecast.rscond.forecast);
% then modify or add some more
Initcond.PAI=PAI;
Initcond.simul_history_end_date=simul_history_end_date;
Initcond.simul_sig=simul_sig;
Initcond.simul_order=simul_order;
Initcond.Qfunc=Qfunc;
Initcond.complementarity=complementarity;
Initcond.simul_shock_uncertainty=obj.options.simul_shock_uncertainty;
Initcond.nsteps=obj.options.simul_periods;
Initcond.k_future=k_future;
Initcond.simul_update_shocks_handle=obj.options.simul_update_shocks_handle;
Initcond.simul_do_update_shocks=obj.options.simul_do_update_shocks;
Initcond.shock_structure=shock_structure;
Initcond.sep_compl=sep_compl;
%-----------------------------------------
if is_conditional_forecasting
    % shocks have already been set from the initial conditions no
    % burn-in simulations necessary
    Initcond.burn=0;
else
    Initcond.burn=obj.options.simul_burn;
    which_shocks=true(1,exo_nbr);
    which_shocks(obj.exogenous.is_observed)=false;
    shocks=utils.forecast.create_shocks(exo_nbr,[],~which_shocks,Initcond);
    y0.econd.data=shocks(:,:,ones(1,3));
end
ncond_regimes=Initcond.nsteps+Initcond.burn;
rcond_data=nan(1,ncond_regimes,3);
if isempty(y0.rcond.data)
    y0.rcond.data=rcond_data;
else
    cutoff=min(ncond_regimes,size(y0.rcond.data,2));
    rcond_data(:,1:cutoff,:)=y0.rcond.data(:,1:cutoff,:);
    y0.rcond.data=rcond_data;
end
simul_regime=obj.options.simul_regime;
if all(isnan(y0.rcond.data(:))) && ~isempty(simul_regime)
    nregs=numel(simul_regime);
    if nregs==1
        y0.rcond.data(:)=simul_regime;
    else
        nregs=min(nregs,size(y0.rcond.data,2));
        y0.rcond.data(1,1:nregs,1)=simul_regime(1:nregs);
        y0.rcond.data=y0.rcond.data(:,:,ones(1,3));
    end
end
Initcond.y=y0;

    function set_shocks_and_states()
        conditional_shocks=obj.options.forecast_cond_exo_vars;
        if ~is_conditional_forecasting
            is_conditional_forecasting=~isempty(conditional_shocks);
        end
        if is_conditional_forecasting
            shock_pos=locate_variables(conditional_shocks,obj.exogenous.name);
        else
            % conditional on all shocks for stochastic simulation, etc.
            shock_pos=1:sum(obj.exogenous.number);
        end
        % shocks
        %---------
        datae=build_cond_data(conditional_shocks);
        y0(1).econd=struct('data',datae,'pos',shock_pos);
        
        % regimes: silence the error
        %----------------------------
        datar=build_cond_data('regime',true);
        y0(1).rcond=struct('data',datar,'pos',nan);
    end

    function set_endogenous_variables()
        
        if has_data
            if do_dsge_var
                endo_names=obj.observables.name(:);
            else
                endo_names=obj.endogenous.name(:);
            end
            endo_names=endo_names(:,ones(1,nlags));
            for ilag=2:nlags
                endo_names(:,ilag)=strcat(endo_names(:,ilag),sprintf('{-%0.0f}',ilag-1));
            end
            endo_names=endo_names(:).';
            
            y0(1).y=utils.forecast.load_start_values(endo_names,...
                simul_historical_data,simul_history_end_date,y0(1).y);
            % keep first page only: conditions are separate
            %-----------------------------------------------
            y0(1).y=y0(1).y(:,:,1);
            
            % past regimes for the computation of initial regime probabilities
            %-----------------------------------------------------------------
            PAI_lag=utils.forecast.load_start_values(obj.markov_chains.regime_names,...
                simul_historical_data,simul_history_end_date,nan(size(PAI)));
            if all(~isnan(PAI_lag))
                PAI=transpose(obj.solution.transition_matrices.Q)*PAI_lag;
            end
        end

        conditional_vars=obj.options.forecast_cond_endo_vars;
        datay=build_cond_data(conditional_vars);
        y_pos=strcmp(conditional_vars,obj.observables.name);
        y_pos=obj.observables.state_id(y_pos);
        y0(1).ycond=struct('data',datay,'pos',y_pos);
        is_conditional_forecasting=~isempty(datay);
    end

    function datax=build_cond_data(names,silent)
        if nargin<2
            silent=false;
        end
        if isempty(names)
            names={};
        elseif ischar(names)
            names=cellstr(names);
        end
        nx=numel(names);
        if has_data
            nspan=simul_historical_data.NumberOfPages-1;
            datax=nan(nx,nspan,3);
            for ivar=1:nx
                vname=names{ivar};
                if silent && ~any(strcmp(vname,simul_historical_data.varnames));
                    di=nan(1,nspan);
                else
                    di=get_data(vname);
                end
                % always hard if not stated otherwise
                %------------------------------------
                datax(ivar,:,1)=di;
                datax(ivar,:,2)=di;
                datax(ivar,:,3)=di;
                % and possibly soft if possible
                %-------------------------------
                apply_soft(['lower_',vname],2)
                apply_soft(['upper_',vname],3)
                % should the user put -inf or inf themselves? most def but I
                % don't believe that if they put a number they will put
                % anything else than finite...
                % - Another interesting case if what happens if the user only
                % puts either the lower bound or just the upper bound.
                % - According to this specification, there is no conditionning
                % on a variable if the user does not give a central tendency.
                % - We have to revisit the role of nan central tendency vs nan
                % bounds...
            end
        else
            nspan=0;
            datax=nan(nx,nspan,3);
        end
        function apply_soft(vname,pos)
            low=any(strcmp(vname,simul_historical_data.varnames));
            if low
                di_=get_data(vname);
                datax(ivar,:,pos)=di_;
            end
        end
        function di=get_data(vname)
            di=simul_historical_data(simul_historical_data.finish,vname);
            di=squeeze(double(di));
            di=di(2:end);
        end
    end
end