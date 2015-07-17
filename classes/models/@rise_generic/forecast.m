function cond_fkst_db=forecast(obj,varargin)
% forecast - computes forecasts for rise|dsge|svar|rfvar models
%
% Syntax
% -------
% ::
%
%   cond_fkst_db=forecast(obj,varargin)
%
% Inputs
% -------
%
% - **obj** [rise|dsge|svar|rfvar]: model object
%
% - **varargin** : additional inputs coming in pairs. These include but are
%   not restricted to:
%
%   - **forecast_to_time_series** [{true}|false]: sets the output to time
%       series format or not
%
%   - **forecast_nsteps** [integer|{12}]: number of forecasting steps
%
%   - **forecast_start_date** [char|numeric|serial date]: date when the
%       forecasts start (end of history + 1)
%
%   - **forecast_conditional_hypothesis** [{jma}|ncp|nas]: in dsge models in
%       which agents have information beyond the current period, this
%       option determines the number of periods of shocks need to match the
%       restrictions:
%       - Hypothesis **jma** assumes that irrespective of how
%           many periods of conditioning information are remaining, agents
%           always receive information on the same number of shocks.
%       - Hypothesis **ncp** assumes there are as many shocks periods as
%           the number of the number of conditioning periods
%       - Hypothesis **nas** assumes there are as many shocks periods as
%           the number of anticipated steps
%
%   - **forecast_cond_endo_vars** [{''},char|cellstr]: names of conditional
%   endogenous variables to be used either in forecasting or in estimation
%
%   - **forecast_cond_exo_vars** [{''},char|cellstr]: names of conditional
%   exogenous variables to be used either in forecasting or in estimation
%
%   - **forecast_shock_uncertainty** [true|{false}]: draw shocks over the
%   simulation horizon.
%
% Outputs
% --------
%
% - **cond_fkst_db** [struct|matrix]: depending on the value of
%   **forecast_to_time_series** the returned output is a structure with
%   time series or a cell containing a matrix and the information to
%   reconstruct the time series.
%
% More About
% ------------
%
% - the historical information as well as the conditioning information come
%   from the same database. The time series must be organized such that for
%   each series, the first page represents the actual data and all
%   subsequent pages represent conditional information. If a particular
%   condition is "nan", that location is not constrained
% - Conditional forecasting for nonlinear models is also supported.
%   However, the solving of the implied nonlinear problem may fail if the
%   model displays instability
% - Both HARD CONDITIONS and SOFT CONDITIONS are implemented but the latter
%   are currently disabled in expectation of a better user interface.
% - The data may also contain time series for a variable with name
%   **regime** in that case, the forecast/simulation paths are computed
%   following the information therein. **regime** must be a member of 1:h,
%   where h is the maximum number of regimes.
%
% Examples
% ---------
%
% See also: simulate

if isempty(obj)
    if nargout>1
        error([mfilename,':: number of output arguments cannot exceed 1 when the object is empty'])
    end
    cond_fkst_db=struct('forecast_conditional_hypothesis',0,...
        'forecast_start_date','',...
        'forecast_shock_uncertainty',false,...
        'forecast_nsteps',12,...
        'forecast_cond_endo_vars','',...
        'forecast_cond_exo_vars','',...
        'forecast_to_time_series',true);
    return
end
%% pass the options
obj=set(obj,varargin{:});

%% pass some options to the simulate method
if isempty(obj.options.forecast_start_date)
    if isempty(obj.options.estim_end_date)
        error('cannot find a date to start the forecast')
    end
    obj=set(obj,'forecast_start_date',serial2date(date2serial(obj.options.estim_end_date)+1));
end

end_date='';
if ~isempty(obj.options.forecast_start_date)
    end_date=serial2date(date2serial(obj.options.forecast_start_date)-1);
end
fkstdata='';
if ~isempty(obj.options.data)
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