%--- help for generic/simulate ---
%
%  simulate - simulates a RISE model
% 
%  ::
% 
% 
%    [db,states,retcode] = simulate(obj,varargin)
% 
%  Args:
% 
%     - **obj** [rfvar|dsge|rise|svar]: model object
% 
%     - **varargin** : additional arguments including but not restricted to
% 
%       - **simul_periods** [integer|{100}]: number of simulation periods
% 
%       - **simul_burn** [integer|{100}]: number of burn-in periods. This
%         should not be confused with forecast_conditional_sampling_burnin,
%         which is used in the sampling from the truncated multivariate
%         normal distribution.
% 
%       - **simul_historical_data** [ts|simplan|struct|{''}]: historical 
%           data from which the simulations are based. If empty, the 
%           simulations start at the steady state.
% 
%       - **simul_history_end_date** [char|integer|serial date]: last date of
%           history
% 
%       - **simul_regime** [integer|vector|function handle|{[]}]: regimes for
%           which the model is simulated. When it is a function handle, then it
%           is a switch_rule and should accept as inputs y (array over all regimes),
%           regimes_1_t_1(regimes from 1 to t-1), sims_1_t_1(the simulated
%           series up to t-1),varargin (possibly further arguments to the
%           function handle). The output is a logical vector that is true for
%           the columns that are acceptable/feasible and false otherwise.
% 
%       - **simul_to_time_series** [{true}|false]: if true, the output is a
%           time series, else a cell array with a matrix and information on
%           elements that help reconstruct the time series.
% 
%       - **simul_shock_uncertainty** [{true}|false]: draw shocks over the
%       simulation horizon.
% 
%       - **simul_anticipate_zero** [true|{false}]: When shocks are drawn, this
%       option allows to impose that agents continue to see only the
%       contemporaneous shocks.
% 
%  Returns: 
% 
%     - **db** [struct|cell array]: if **simul_to_time_series** is true, the
%       output is a time series, else a cell array with a matrix and
%       information on elements that help reconstruct the time series.
% 
%     - **states** [vector]: history of the regimes over the forecast horizon
% 
%     - **retcode** [integer]: if 0, the simulation went fine. Else something
%       got wrong. In that case one can understand the problem by running
%       decipher(retcode)
% 
%  Note:
% 
%     - **simul_historical_data** contains the historical data as well as
%       conditional information over the forecast horizon. It may also include
%       as an alternative to **simul_regime**, a time series with name
%       **regime**, which indicates the regimes over the forecast horizon.
% 
%  Example:
% 
%     See also:
%
%    Other uses of simulate
%
%       arima/simulate             lassoblm/simulate
%       bates/simulate             merton/simulate
%       bssm/simulate              mixconjugateblm/simulate
%       cir/simulate               mixsemiconjugateblm/simulate
%       conjugateblm/simulate      msVAR/simulate
%       conjugatebvarm/simulate    normalbvarm/simulate
%       customblm/simulate         quantumCircuit/simulate
%       diffuseblm/simulate        regARIMA/simulate
%       diffusebvarm/simulate      sde/simulate
%       dsge/simulate              semiconjugateblm/simulate
%       dtmc/simulate              semiconjugatebvarm/simulate
%       egarch/simulate            ssm/simulate
%       empiricalblm/simulate      tsVAR/simulate
%       garch/simulate             varm/simulate
%       gjr/simulate               vecm/simulate
%       heston/simulate
%