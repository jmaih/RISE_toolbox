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
%       should not be confused with forecast_conditional_sampling_burnin, which
%       is used in the sampling from the truncated multivariate normal
%       distribution.
% 
%       - **simul_historical_data** [ts|struct|{''}]: historical data from
%           which the simulations are based. If empty, the simulations start at
%           the steady state.
% 
%       - **simul_history_end_date** [char|integer|serial date]: last date of
%           history
% 
%       - **simul_regime** [integer|vector|function handle|{[]}]: regimes for
%           which the model is simulated. When it is a function handle, then it
%           should accept as inputs y (array over all regimes),
%           regimes_1_t_1(regimes from 1 to t-1), sims_1_t_1(the simulated
%           series up to t-1),varargin (possibly further arguments to the
%           function handle). The output is a logical vector that is true for
%           the columns that are acceptable/feasible and false otherwise.
% 
%       - **simul_to_time_series** [{true}|false]: if true, the output is a
%           time series, else a cell array with a matrix and information on
%           elements that help reconstruct the time series.
% 
%       - **simul_honor_constraints** [true|{false}]: honor restrictions during
%           simulations. If true, agents have to be able anticipate the future.
% 
%       - **simul_frwrd_back_shoot** [true|{false}]: uses an algorithm that
%       checks the constraints are satisfied at all horizons instead of just
%       one at a time.
% 
%       - **simul_shock_uncertainty** [{true}|false]: draw shocks over the
%       simulation horizon.
% 
%       - **simul_honor_constraints_through_switch** [true|{false}]: if true,
%       constraints are honored through the switching mechanism. In that case
%       the number of regimes should be greater than 1. If false, constraints
%       are honored through an anticipatory behavior. In this case, there
%       should be shocks that are foreseen.
% 
%       - **simul_anticipate_zero** [true|{false}]: When shocks are drawn, this
%       option allows to impose that agents continue to see only the
%       contemporaneous shocks.
% 
%       - **simul_bgp_deviation** [true|{false}]: When the model is
%       nonstationary, a growth component appears in the solution. This option
%       enables or disables that component.
% 
%  Returns: 
%     :
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
%    Other functions named simulate
%
%       arima/simulate           garch/simulate
%       conjugateblm/simulate    gjr/simulate
%       customblm/simulate       regARIMA/simulate
%       diffuseblm/simulate      sde/simulate
%       dsge/simulate            semiconjugateblm/simulate
%       dtmc/simulate            ssm/simulate
%       egarch/simulate          varm/simulate
%       empiricalblm/simulate    vecm/simulate
%