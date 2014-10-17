function Initcond=set_simulation_initial_conditions(obj)
% H1 line
%
% Syntax
% -------
% ::
%
% Inputs
% -------
%
% Outputs
% --------
%
% More About
% ------------
%
% Examples
% ---------
%
% See also: 



% complementarity: do it here so as to keep the memory light
%-----------------------------------------------------------
if ~isfield(obj.routines,'complementarity')||... I do not expect VARs to have this although in theory they could.
        isempty(obj.routines.complementarity)
    complementarity=@(varargin)true;
else
    complementarity=@build_complementarity;
end

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
    'simul_do_update_shocks',obj.options.simul_do_update_shocks);
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

%-----------------------------------------
Initcond.states=states;
Initcond.shocks=shocks;

    function c=build_complementarity(y)
        c=true;
        param=obj.parameter_values(:,1);
        x=[];
        ss=obj.solution.ss{1};
        sparam=[];
        def=obj.solution.definitions{1};
        s0=1;
        s1=1;
        for ii=1:numel(obj.routines.complementarity)
            c= c && obj.routines.complementarity{ii}(y,x,ss,param,sparam,def,s0,s1); 
            if ~c
                break
            end
        end
    end

    function set_shocks_and_states()
        % check there is no shock with name regime (this should be done right
        % from the parser)
        %--------------------------------------------------------------------
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
                if any(isnan(shocks(end,:)))
                    warning('the nan locations in the regimes will be chosen randomly according to the transition probabilities')
                end
%                 for icol=1:nperiods+k
%                     if isnan(shocks(end,icol))
%                         if icol==1
%                             % set the first period to the first regime
%                             shocks(end,icol)=1;
%                             warning('the first regime was nan and has been set to 1 in the simulations')
%                         else
%                             shocks(end,icol)=shocks(end,icol-1);
%                         end
%                     end
%                 end
                % make sure there is no regime exceeding h and all are positive integers
                %-----------------------------------------------------------------------
                regimes_row=shocks(end,:);
                shocks(end,:)=[];
                h=obj.markov_chains.regimes_number;
                good=~isnan(regimes_row);
                if ~(all(ceil(regimes_row(good))==floor(regimes_row(good))) && ...
                        all(regimes_row(good)>=1) &&...
                        all(regimes_row(good)<=h))
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
        endo_names=obj.endogenous.name(:);
        endo_names=endo_names(:,ones(1,nlags));
        for ilag=2:nlags
            endo_names(:,ilag)=strcat(endo_names(:,ilag),sprintf('{-%0.0f}',ilag-1));
        end
        endo_names=endo_names(:).';

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