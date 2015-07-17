function myirfs=irf(obj,varargin)
% irf - computes impulse responses for a RISE model
%
% Syntax
% -------
% ::
%
%   myirfs=irf(obj)
%
%   myirfs=irf(obj,varargin)
%
% Inputs
% -------
%
% - **obj** [rise|dsge|rfvar|svar]: single or vector of RISE models
%
% - **varargin** : optional options coming in pairs. The notable ones that
%   will influence the behavior of the impulse responses are:
%
% - **irf_shock_list** [char|cellstr|{''}]: list of shocks for which we
%   want to compute impulse responses
%
% - **irf_var_list** [char|cellstr|{''}]: list of the endogenous variables
%   we want to report
%
% - **irf_periods** [integer|{40}]: length of the irfs
%
% - **irf_shock_sign** [numeric|-1|{1}]: sign or scale of the original
%   impulse. If **irf_shock_sign** >0, we get impulse responses to a
%   positive shock. If **irf_shock_sign** <0, the responses are negative.
%   If If **irf_shock_sign** =0, all the responses are 0.
%
% - **irf_draws** [integer|{50}]: number of draws used in the simulation
%   impulse responses in a nonlinear model. A nonlinear model is defined as
%   a model that satisfies at least one of the following criteria
%   - solved at an order >1
%   - has more than one regime and option **irf_regime_specific** below is
%       set to false
%
% - **irf_type** [{irf}|girf]: type of irfs. If the type is irf, the
%   impulse responses are computed directly exploiting the fact that the
%   model is linear. If the type is girf, the formula for the generalized
%   impulse responses is used: the irf is defined as the expectation of the
%   difference of two simulation paths. In the first path the initial
%   impulse for the shock of interest is nonzero while it is zero for the
%   second path. All other shocks are the same for both paths in a given
%   simulation.
%
% - **irf_regime_specific** [{true}|false]: In a switching model, we may or
%   may not want to compute impulse responses specific to each regime.
%
% - **irf_use_historical_data** [{false}|true]: if true, the data stored in
%   option **simul_historical_data** are used as initial conditions. But
%   the model has to be nonlinear otherwise the initial conditions are set
%   to zero. This option gives the flexibility to set the initial
%   conditions for the impulse responses.
%
% - **irf_to_time_series** [{true}|false]: If true, the output is in the
%   form of time series. Else it is in the form of a cell containing the
%   information needed to reconstruct the time series.
%
% Outputs
% --------
%
% - **myirfs** [{struct}|cell]: Impulse response data
%
% More About
% ------------
%
% - for linear models or models solved up to first order, the initial
%   conditions as well as the steady states are set to 0 in the computation
%   of the impulse responses.
%
% - for nonlinear models, the initial conditions is the ergodic mean
%
% Examples
% ---------
%
% See also: 

too_small=1e-9;

num_fin=@(x)isnumeric(x) && isscalar(x) && isfinite(x);
num_fin_int=@(x)num_fin(x) && floor(x)==ceil(x) && x>=0;
mydefaults={'irf_shock_list','',@(x)ischar(x)||iscellstr(x),'irf_shock_list must be char or cellstr'
        'irf_var_list','',@(x)ischar(x)||iscellstr(x),'irf_var_list must be char or cellstr'
        'irf_periods',40,@(x)num_fin_int(x),'irf_periods must be a finite and positive integer'
        'irf_shock_sign',1,@(x)num_fin(x) && x~=0,'irf_shock_sign must be a numeric scalar different from 0'
        'irf_draws',50,@(x)num_fin_int(x),'irf_draws must be a finite and positive integer'
        'irf_type','irf',@(x)any(strcmp(x,{'irf','girf'})),'irf_type must be irf or girf'
        'irf_regime_specific',true,@(x)islogical(x),'irf_regime_specific must be a logical'
        'irf_use_historical_data',false,@(x)islogical(x),'irf_use_historical_data must be a logical'
        'irf_to_time_series',true,@(x)islogical(x),'irf_to_time_series must be a logical'
        };
if isempty(obj)
    myirfs=cell2struct(mydefaults(:,2),mydefaults(:,1),1);
    return
end

nobj=numel(obj);
obj=set(obj,varargin{:});

myirfs=cell(1,nobj);
% check that the models are consistent
check_irf_consistency(obj)
for ii=1:nobj
    myirfs{ii}=irf_intern(obj(ii));
end

myirfs=format_irf_output(myirfs);

    function dsge_irfs=irf_intern(obj)
        if nobj>1
            obj.options.irf_to_time_series=true;
        end
        [irf_shock_list,irf_var_list,irf_periods,irf_shock_sign,irf_draws,...
            irf_type,irf_regime_specific,irf_use_historical_data,...
            irf_to_time_series]=utils.miscellaneous.parse_arguments(mydefaults,...
            'irf_shock_list',obj.options.irf_shock_list,...
            'irf_var_list',obj.options.irf_var_list,...
            'irf_periods',obj.options.irf_periods,...
            'irf_shock_sign',obj.options.irf_shock_sign,...
            'irf_draws',obj.options.irf_draws,...
            'irf_type',obj.options.irf_type,...
            'irf_regime_specific',obj.options.irf_regime_specific,...
            'irf_use_historical_data',obj.options.irf_use_historical_data,...
            'irf_to_time_series',obj.options.irf_to_time_series);
        is_dsge=isa(obj,'dsge');
        
        obj.options.simul_periods=irf_periods;
        obj.options.simul_burn=0;
        if ~irf_use_historical_data
            obj.options.simul_historical_data=ts.empty(0);
            obj.options.simul_history_end_date='';
        end
        exo_nbr=sum(obj.exogenous.number);
        
        if isempty(irf_var_list)
            if is_dsge
                irf_var_list=get(obj,'endo_list(original)');
            else
                irf_var_list=get(obj,'endo_list');
            end
        elseif ischar(irf_var_list)
            irf_var_list=cellstr(irf_var_list);
        end
        
        detList=get(obj,'exo_list(observed)');
        exoList=get(obj,'exo_list');
        if isempty(irf_shock_list)
            irf_shock_list=get(obj,'exo_list(~observed)');
        end
        if ischar(irf_shock_list)
            irf_shock_list=cellstr(irf_shock_list);
        end
        % sort in case the user has not listed the shocks in the alphabetic
        % order
        %------------------------------------------------------------------
        irf_shock_list=sort(irf_shock_list);
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
        if any(ismember(irf_shock_list,get(obj,'exo_list(observed)')))
            error('cannot compute irfs of observed shocks')
        end
        which_shocks=false(1,exo_nbr);
        which_shocks(position)=true;
        nshocks=sum(which_shocks);
        det_shocks=false(1,exo_nbr);
        det_pos=locate_variables(detList,exoList,true);
        det_shocks(det_pos)=true;
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
        do_dsge_var=false;
        if is_dsge
            do_dsge_var=obj.is_dsge_var_model && obj.options.dsgevar_var_regime;
            solve_order=obj.options.solve_order;
            % hide future shocks if required
            %-------------------------------
            obj=do_not_anticipate_future_shocks(obj);
        end
        % load the order_var solution
        %-----------------------------
        [T,~,steady_state,new_order,state_vars_location]=load_solution(obj,'ov',do_dsge_var);
       
        % initial conditions
        %-------------------
        Initcond=set_simulation_initial_conditions(obj);
        
        h=obj.markov_chains.regimes_number;
        
        Initcond=utils.forecast.initial_conditions_to_order_var(Initcond,new_order,obj.options);
        
        if solve_order==1 % first-order solution ||isa(obj,'rfvar')
            % kill the steady state and the initial conditions
            for ireg=1:h
                steady_state{ireg}=0*steady_state{ireg};
            end
            Initcond.y.y=0*Initcond.y.y;
        end
        y0=Initcond.y;
                   
        girf=solve_order>1||(solve_order==1 && h>1 && strcmp(irf_type,'girf'));
        if ~girf
            irf_draws=1;
        end
        if ~obj.options.irf_use_historical_data
            Initcond.shocks=0*Initcond.shocks;
        end
        irf_shock_uncertainty=irf_draws>1;
        number_of_threads=h;
        if ~irf_regime_specific
            number_of_threads=1;
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
        
        % initialize output: use imaginaries for the non-computed
        % variables. This is so that the time series object does not squeak
        %------------------------------------------------------------------
        endo_nbr=obj.endogenous.number;
        Impulse_dsge=zeros(endo_nbr,Initcond.nsteps+1,nshocks,irf_draws,...
            number_of_threads)+1i;
        % select only the relevant rows in case we are dealing with
        % a VAR with many lags a BVAR_DSGE
        %----------------------------------------------------------
        relevant=1:endo_nbr;
        max_rows=endo_nbr;
        if do_dsge_var
            max_rows=obj.observables.number(1);
            relevant=obj.inv_order_var(obj.observables.state_id);
        end
        retcode=0;
        for istate=1:number_of_threads
            if ~retcode
                if h==1||number_of_threads==h
                    Initcond.states(:,1)=istate;
                end
                [xxxx,retcode]=utils.forecast.irf(y0,T,steady_state,...
                    state_vars_location,which_shocks,det_shocks,Initcond);
                Impulse_dsge(relevant,:,:,:,istate)=xxxx(1:max_rows,:,:,:);
            end
        end
        
        % set to 0 the terms that are too tiny
        Impulse_dsge(abs(Impulse_dsge)<=too_small)=0;
        
        % re-order the variables according to the inv_order_var;
        Impulse_dsge=re_order_output_rows(obj,Impulse_dsge);
        % going from variables x time x shocks x irf_draws x regimes
        % reshape as time x regimes x variables x shocks x irf_draws
        %------------------------------------------------------------
        Impulse_dsge=permute(Impulse_dsge,[2,5,1,3,4]);
        
        % average across  the irf_draws dimension
        %----------------------------------------
        Impulse_dsge=mean(Impulse_dsge,5);
        
        % distribution of irfs
        %---------------------
        dsge_irfs=distribute_irfs();
        
        function dsge_irfs=distribute_irfs()
            startdate=0;
            if number_of_threads>1
                RegimeNames=strcat('regime_',num2str((1:h)'));
            else
                RegimeNames=irf_type;
            end
            RegimeNames=cellfun(@(x)x(~isspace(x)),cellstr(RegimeNames),'uniformOutput',false);
            if irf_to_time_series
                dsge_irfs=struct();
                vlocs=locate_variables(irf_var_list,get(obj,'endo_list'));
                for ishock=1:nshocks
                    shock_name=irf_shock_list{ishock};
                    for vv=1:numel(irf_var_list)
                        dsge_irfs.(shock_name).(irf_var_list{vv})=...
                            ts(startdate,squeeze(Impulse_dsge(:,:,vlocs(vv),ishock)),RegimeNames);
                    end
                end
            else
                dsge_irfs={Impulse_dsge,...
                    {
                    '1=time',Initcond.nsteps+1
                    '2=regime names',RegimeNames
                    '3= endogenous names',irf_var_list
                    '4= shock names',irf_shock_list
                    }
                    };
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
