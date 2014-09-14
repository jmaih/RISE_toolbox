function Initcond=set_simulation_initial_conditions(obj)

regimes_number=obj.markov_chains.regimes_number;
PAI=1/regimes_number*ones(regimes_number,1); % ErgodicDistribution(Q)

y0=struct('y',{},'y_lin',{});
nlags=1;
if isa(obj,'svar')
    nlags=obj.nlags;
elseif isa(obj,'dsge')
else
    error(['model of class ',class(obj),' not ready for simulation'])
end

for ireg=1:regimes_number
    y0(ireg).y=vec(obj.solution.ss{ireg}(:,ones(1,nlags)));
end

simul_history_end_date=0;
simul_historical_data=obj.options.simul_historical_data;
shocks=[];
states=[];
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
end

Q={obj.solution.transition_matrices.Q,[],[]};
if isa(obj,'dsge') && obj.is_endogenous_switching_model
    % take a handle to a private function. we can't access it otherwise
    Q{2}=obj.routines.transition_matrix;
    M=obj.parameter_values;
    % transition matrix should be invariant. And so, hopefully, the first
    % argument to Q{2} could be the updated data in any state.
    % order of the input arguments is y,x,ss,param,sparam,def,s0,s1 
    % remaining arguments of Q{2} after the first one
    %-------------------------------------------------
    Q{3}={[],obj.solution.ss{1},mean(M,2),[],[],[],[]};
end
simul_pruned=false;
simul_sig=0;
simul_order=1;
k_future=0;
if isa(obj,'dsge')
    simul_sig=obj.options.simul_sig;
    simul_pruned=obj.options.simul_pruned;
    if isempty(obj.options.simul_order);
        simul_order=obj.options.solve_order;
    else
        simul_order=obj.options.simul_order;
    end
    k_future=max(obj.exogenous.shock_horizon);
end
if ~simul_pruned
    y0=rmfield(y0,'y_lin');
end
Initcond=struct('y',{y0},...
    'PAI',PAI,...
    'simul_history_end_date',simul_history_end_date,...
    'simul_sig',simul_sig,...
    'simul_order',simul_order,...
    'Q',{Q},...
    'random',true,...
    'nsteps',obj.options.simul_periods,...
    'k_future',k_future);
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
    shocks=utils.forecast.create_shocks(exo_nbr,[],~which_shocks,Initcond);
end
if isempty(states)
    states=nan(Initcond.nsteps+Initcond.burn,1);
end

%-----------------------------------------
Initcond.states=states;
Initcond.shocks=shocks;
Initcond.states=states;

    function set_shocks_and_states()
        % check there is no shock with name regime (this should be done right
        % from the parser)
        %--------------------------------------------------------------------
        warning('the lines below should be moved to the parser')
        if any(strcmp(obj.exogenous.name,'regime'))
            error('no exogenous variable can be called regime')
        end
        new_shock_names=[obj.exogenous.name,'regime'];
        
        varnames=simul_historical_data.varnames;
        if any(ismember(new_shock_names,varnames))
            % locate the left_date and right_date in the historical database dates
            %---------------------------------------------------------------------
            [right,flag_r]=date2obs(simul_historical_data.start,...
                date2serial(simul_history_end_date)+1);
            if flag_r
                warning(['too few observations to define initial conditions for shocks. ',...
                    'Maybe you should adjust simul_history_end_date?'])
            else
                raw_data=double(simul_historical_data).';
                % This exo_nbr includes the place holder for the regime
                %------------------------------------------------------
                exo_plus_regime_nbr=numel(new_shock_names);
                k=0;
                if isa(obj,'dsge')
                    k=max(obj.exogenous.shock_horizon);
                end
                nperiods=size(raw_data(:,right:end),2);
                % initialize shocks and nan the regime row
                %-----------------------------------------
                simul_periods=obj.options.simul_periods;
                shocks=zeros(exo_plus_regime_nbr,max(simul_periods,nperiods)+k);
                shocks(end,:)=nan;
                exo_locs=locate_variables(new_shock_names,varnames,true);
                good=~isnan(exo_locs);
                shocks(good,1:nperiods)=raw_data(exo_locs(good),right:end);
                % make sure the regime row does not have nans
                %--------------------------------------------
                for icol=1:nperiods+k
                    if isnan(shocks(end,icol))
                        if icol==1
                            % set the first period to the first regime
                            shocks(end,icol)=1;
                            warning('the first regime was nan and has been set to 1 in the simulations')
                        else
                            shocks(end,icol)=shocks(end,icol-1);
                        end
                    end
                end
                % make sure there is no regime exceeding h and all are positive integers
                %-----------------------------------------------------------------------
                regimes_row=shocks(end,:);
                shocks(end,:)=[];
                h=obj.markov_chains.regimes_number;
                if ~all(ismember(regimes_row,1:h))
                    error(['regimes must be positive integers and cannot exceed ',int2str(h)])
                end
                % zero all nans in the shocks. In a conditional forecasting exercise,
                % the nan locations have to be found, but this is not what we are doing
                % here under simulation
                shocks(isnan(shocks))=0;
                
                % Now load the states
                %--------------------
                if any(regimes_row)
                    states=regimes_row(:);
                end
                
                % Now load the regimes probabilities
                %-----------------------------------
                
            end
        end
    end

    function set_endogenous_variables()
        endo_names=obj.endogenous.name;
        tmp=endo_names;
        for ilag=1:nlags-1
            endo_names=[endo_names,strcat(tmp,sprintf('{-%0.0f}',ilag))];
        end
        %--------------------------------
        y00=y0(1).y;
        y0(1).y=utils.forecast.load_start_values(endo_names,simul_historical_data,...
            simul_history_end_date,y0(1).y);
        is_changed=abs(y00-y0(1).y)>sqrt(eps);
        for istate=2:regimes_number
            % replace only the changed locations
            y0(istate).y(is_changed)=y0(1).y(is_changed);
        end
    end
end