%--- help for dsge/bvar_dsge ---
%
%  INTERNAL FUNCTION: intermediary file for computing key elements for dsge-var
% 
%  ::
% 
%    [obj,retcode] = bvar_dsge(obj,varargin)
% 
%  Args:
% 
%     obj (rise | dsge): model object
%     varargin : ususal dsge options. The most important of which are:
% 
%       - **dsgevar_lag** [integer|{4}]: number of lags in the VAR
%       - **dsgevar_constant** [false|{true}]: flag for having a constant in
%         the VAR
%       - **dsgevar_var_regime** [false|{true}]: use the VAR in simulations.
%         Otherwise use the DSGE
%       - **dsgevar_inner_param_uncertainty** [true|{false}]: make random draws
%         for the parameters around the mode for each simulation
% 
%  Returns:
%     :
% 
%     - **obj** [rise|dsge]: model object
% 
%  Note:
% 
%     - Because the BVAR-DSGE will have fewer variables than the DSGE, in
%       simulations, the missing variables will be stored as 0+1i in time series.
%       This is true for forecast, irf and simulate
% 
%