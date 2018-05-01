function cond_fkst_db=forecast(obj,varargin)
% forecast - computes forecasts for rise|dsge|svar|rfvar models
%
% ::
%
%
%   cond_fkst_db=forecast(obj,varargin)
%
% Args:
%
%    - **obj** [rise|dsge|svar|rfvar]: model object
%
%    - **varargin** : additional inputs coming in pairs. These include but are
%      not restricted to:
%
%      - **forecast_to_time_series** [{true}|false]: sets the output to time
%          series format or not
%
%      - **forecast_nsteps** [integer|{12}]: number of forecasting steps
%
%      - **forecast_start_date** [char|numeric|serial date]: date when the
%          forecasts start (end of history + 1)
%
%      - **forecast_cond_endo_vars** [{''},char|cellstr]: names of conditional
%      endogenous variables to be used either in forecasting or in estimation
%
%      - **forecast_cond_exo_vars** [{''},char|cellstr]: names of conditional
%      exogenous variables to be used either in forecasting or in estimation
%
%      - **forecast_endo_exo_vars** [{''},char|cellstr]: names of exogenous
%      variables that are potentially modified under conditional forecasting
%      with forward-back shooting
%
%      - **forecast_shock_uncertainty** [true|{false}]: draw shocks over the
%      simulation horizon.
%
% Returns:
%    :
%
%    - **cond_fkst_db** [struct|matrix]: depending on the value of
%      **forecast_to_time_series** the returned output is a structure with
%      time series or a cell containing a matrix and the information to
%      reconstruct the time series.
%
% Note:
%
%    - the historical information as well as the conditioning information come
%      from the same database. The time series must be organized such that for
%      each series, the first page represents the actual data and all
%      subsequent pages represent conditional information. If a particular
%      condition is "nan", that location is not constrained
%
%    - Conditional forecasting for nonlinear models is also supported.
%      However, the solving of the implied nonlinear problem may fail if the
%      model displays instability
%
%    - Both HARD CONDITIONS and SOFT CONDITIONS are implemented. In order to
%    do soft conditions, the variables one wishes to condition on must have
%    lower and upper bounds represented by the presence of variables with
%    names "lower_CONDNAME" and "upper_CONDNAME", where CONDNAME is the name
%    of a particular variable we want to condition on.
%
%    - The data may also contain time series for a variable with name
%      **regime** in that case, the forecast/simulation paths are computed
%      following the information therein. **regime** must be a member of 1:h,
%      where h is the maximum number of regimes.
%
%    - Further options for conditional forecasting can be found in
%    RSCF.FORECAST
%
% Example:
%
%    See also: RISE_GENERIC/SIMULATE, RSCF.FORECAST

if isempty(obj)
    
    mydefaults=the_defaults();
    
    mydefaults=[mydefaults
        utils.forecast.rscond.forecast()];
    
    % the following options are set elsewhere
    badrows=ismember(mydefaults(:,1),{'debug','simul_shock_uncertainty'});
    
    mydefaults(badrows,:)=[];
    
    if nargout
        
        cond_fkst_db=mydefaults;
        
    else
        
        disp_defaults(mydefaults);
        
    end
    
    return
    
end

% pass the options
%--------------------
obj=set(obj,varargin{:});

% pass some options to the simulate method
%--------------------------------------------
if isempty(obj.options.forecast_start_date)
    
    if isempty(obj.options.estim_end_date)
        
        error('cannot find a date to start the forecast')
        
    end
    
    obj=set(obj,'forecast_start_date',serial2date(date2serial(obj.options.estim_end_date)+1));

end

end_date=obj.options.simul_history_end_date;

if ~isempty(obj.options.forecast_start_date)
    
    end_date=serial2date(date2serial(obj.options.forecast_start_date)-1);
    
end

fkstdata=obj.options.simul_historical_data;

if isempty(fkstdata) && ~isempty(obj.options.data)
    
    fkstdata=obj.options.data;
    
end

obj=set(obj,...
    'simul_periods',obj.options.forecast_nsteps,...
    'simul_history_end_date',end_date,...
    'simul_historical_data',fkstdata,...
    'simul_burn',0,...
    'simul_shock_uncertainty',obj.options.forecast_shock_uncertainty);

obj.options.simul_to_time_series=obj.options.forecast_to_time_series;

cond_fkst_db=simulate(obj);

end

function d=the_defaults()

num_fin=@(x)isnumeric(x) && isscalar(x) && isfinite(x);

num_fin_int=@(x)num_fin(x) && floor(x)==ceil(x) && x>=0;

d={
    'forecast_cond_endo_vars','',@(x)ischar(x)||iscellstr(x),...
    'forecast_cond_endo_vars must be char or cellstr'
    
    'forecast_cond_exo_vars','',@(x)ischar(x)||iscellstr(x),...
    'forecast_cond_exo_vars must be char or cellstr'
    
    'forecast_endo_exo_vars','',@(x)ischar(x)||iscellstr(x),...
    'forecast_endo_exo_vars must be char or cellstr'
    
    'forecast_to_time_series',true,@(x)islogical(x),...
    'forecast_to_time_series must be a logical'
    
    'forecast_nsteps',12,@(x)num_fin_int(x) && x>0,...
    'forecast_nsteps must be a finite and positive integer'
    
    'forecast_shock_uncertainty',false,@(x)islogical(x),...
    'forecast_shock_uncertainty must be a logical'
    
    'forecast_start_date','',@(x)is_date(x)||is_serial(x),...
    'forecast_start_date must be a valid date'
    };

end