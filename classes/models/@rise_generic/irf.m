function dsge_irfs=irf(obj,varargin)

too_small=1e-9;

if isempty(obj)
    dsge_irfs=struct('irf_shock_list','',...
        'irf_var_list','',...
        'irf_periods',40,...
        'irf_shock_sign',1,...
        'irf_draws',50,...
        'irf_type','irf',...
        'irf_regime_specific',true);
    return
end

nobj=numel(obj);

dsge_irfs=cell(1,nobj);
% check that the models are consistent
check_irf_consistency(obj)
for ii=1:nobj
    dsge_irfs{ii}=irf_intern(obj(ii));
end

dsge_irfs=format_irf_output(dsge_irfs);

    function dsge_irfs=irf_intern(obj)
        
        obj=set(obj,varargin{:});
        is_dsge=isa(obj,'dsge');
        
        obj.options.simul_periods=obj.options.irf_periods	  ;
        obj.options.simul_burn=0;
        obj.options.simul_historical_data=ts.empty(0);
        obj.options.simul_history_end_date='';
        irf_shock_list        =obj.options.irf_shock_list;
        irf_var_list          =obj.options.irf_var_list  ;
        irf_shock_sign        =obj.options.irf_shock_sign;
        irf_type	          =obj.options.irf_type	  ;
        irf_draws	          =obj.options.irf_draws	  ;
        irf_regime_specific   =obj.options.irf_regime_specific;
        exo_nbr=sum(obj.exogenous.number);
        which_shocks=true(1,exo_nbr);
        which_shocks(obj.exogenous.is_observed)=false;
        
        if isempty(irf_var_list)
            if is_dsge
                irf_var_list=get(obj,'endo_list(original)');
            else
                irf_var_list=get(obj,'endo_list');
            end
        elseif ischar(irf_var_list)
            irf_var_list=cellstr(irf_var_list);
        end
        
        exoList=get(obj,'exo_list(~observed)');
        if isempty(irf_shock_list)
            irf_shock_list=exoList;
        end
        if ischar(irf_shock_list)
            irf_shock_list=cellstr(irf_shock_list);
        end
        position=locate_variables(irf_shock_list,exoList,true);
        if any(isnan(position))
            disp(irf_shock_list(isnan(position)))
            if isa(obj,'rfvar')
                error(['The reduced-form var has been identified. ',...
                    'List the structural shocks instead of the reduced-form ones'])
            else
                error('The above list of shocks cannot be used in irf')
            end
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
        solve_order=1;
        if is_dsge
            solve_order=obj.options.solve_order;
            % hide future shocks if required
            %-------------------------------
            obj=do_not_anticipate_future_shocks(obj);
        end
        nshocks=numel(irf_shock_list);
        
        % initial conditions
        %-------------------
        Initcond=generic_tools.set_simulation_initial_conditions(obj);
        
        h=obj.markov_chains.regimes_number;
        % load the order_var solution
        %-----------------------------
        [T,~,steady_state,new_order,state_vars_location]=load_solution(obj,'ov');
        y0=Initcond.y;
        for ireg=1:h
            y0(ireg).y=y0(ireg).y(new_order,:);
        end
        
        girf=solve_order>1||(solve_order==1 && h>1 && strcmp(irf_type,'girf'));
        %         quash_regimes=(h>1 && ~irf_regime_specific);
        %||(solve_order==1 && girf)
        if ~girf
            irf_draws=1;
        end
        
        irf_shock_uncertainty=irf_draws>1;
        number_of_threads=h;
        if ~irf_regime_specific
            number_of_threads=1;
        end
        if number_of_threads==1
            %             [y0,T,steady_state]=...
            %                 utils.forecast.aggregate_initial_conditions(Initcond.PAI,...
            %                 y0,T,steady_state);
            y0=utils.forecast.aggregate_initial_conditions(Initcond.PAI,y0);
        end
        
        % further options
        %----------------
        further_options={
            'nsimul',irf_draws
            'impulse',1*irf_shock_sign
            'random',irf_shock_uncertainty
            'girf',girf
            };
        for irow=1:size(further_options,1)
            opname=further_options{irow,1};
            opval=further_options{irow,2};
            Initcond.(opname)=opval;
        end
        % the shocks drawn in the initial conditions will be ignored
        
        % initialize output
        %------------------
        endo_nbr=obj.endogenous.number(end);
        Impulse_dsge=zeros(endo_nbr,Initcond.nsteps+1,nshocks,irf_draws,number_of_threads);
        retcode=0;
        for istate=1:number_of_threads
            if ~retcode
                if h==1||number_of_threads==h
                    Initcond.states(:,1)=istate;
                end
                [xxxx,retcode]=utils.forecast.irf(y0(istate),T,steady_state,...
                    state_vars_location,which_shocks,Initcond);
                % select only the relevant rows in case we are dealing with
                % a VAR with many lags
                %----------------------------------------------------------
                Impulse_dsge(:,:,:,:,istate)=xxxx(1:endo_nbr,:,:,:);
            end
        end
        
        % set to 0 the terms that are too tiny
        Impulse_dsge(abs(Impulse_dsge)<=too_small)=0;
        
        % re-order the variables according to the inv_order_var;
        if is_dsge
            Impulse_dsge=Impulse_dsge(obj.inv_order_var.after_solve,:,:,:,:);
        end
        % going from variables x time x shocks x irf_draws x regimes
        % reshape as time x regimes x variables x shocks x irf_draws
        %------------------------------------------------------------
        Impulse_dsge=permute(Impulse_dsge,[2,5,1,3,4]);
        
        % average across  the irf_draws dimension
        %----------------------------------------
        Impulse_dsge=mean(Impulse_dsge,5);
        
        % distribution of irfs
        %---------------------
        startdate=0;
        if number_of_threads>1
            RegimeNames=strcat('regime_',num2str((1:h)'));
        else
            RegimeNames=irf_type;
        end
        RegimeNames=cellfun(@(x)x(~isspace(x)),cellstr(RegimeNames),'uniformOutput',false);
        dsge_irfs=struct();
        vlocs=locate_variables(irf_var_list,get(obj,'endo_list'));
        for ishock=1:exo_nbr
            shock_name=irf_shock_list{ishock};
            for vv=1:numel(irf_var_list)
                dsge_irfs.(shock_name).(irf_var_list{vv})=...
                    ts(startdate,squeeze(Impulse_dsge(:,:,vlocs(vv),ishock)),RegimeNames);
            end
        end
        
    end
end

function check_irf_consistency(obj)
nobj=numel(obj);
if nobj>1
    first_list={'endogenous','exogenous'};
    second_list={'regimes'};
    third_list={'irf_shock_list','irf_var_list','irf_periods','irf_type'};
    for ilist=1:numel(first_list)
        ref_list=obj(1).(first_list{ilist}).name;
        for iobj=2:nobj
            list=obj(iobj).(first_list{ilist}).name;
            if ~isequal(ref_list,list)
                warning([first_list{ilist},' is not the same across models'])
            end
        end
    end
    for ilist=1:numel(second_list)
        ref_list=obj(1).markov_chains.(second_list{ilist});
        for iobj=2:nobj
            list=obj(iobj).markov_chains.(second_list{ilist});
            if ~isequal(ref_list,list)
                warning([second_list{ilist},' is not the same across models'])
            end
        end
    end
    for ilist=1:numel(third_list)
        ref_list={obj(1).options.(third_list{ilist})};
        for iobj=2:nobj
            list={obj(iobj).options.(third_list{ilist})};
            if ~isequal(ref_list,list)
                error([third_list{ilist},' should be the same across models'])
            end
        end
    end
end
end

function dsge_irfs=format_irf_output(dsge_irfs)
nobj=numel(dsge_irfs);
if nobj==1
    dsge_irfs=dsge_irfs{1};
else
    if isempty(dsge_irfs{1})
        return
    end
    shockList=fieldnames(dsge_irfs{1});
    tmp=struct();
    shock_models=cell(1,nobj);
    for ishock=1:numel(shockList)
        for mm=1:nobj
            shock_models{mm}=dsge_irfs{mm}.(shockList{ishock});
        end
        tmp.(shockList{ishock})=utils.time_series.concatenate_series_from_different_models(shock_models);
    end
    % aggregate
    dsge_irfs=tmp;
end
end
