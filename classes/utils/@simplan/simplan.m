%  `simplan` - Simulation and forecasting plan for DSGE models.
% 
%  Syntax:
%    obj = simplan(model, horizon, initialRegime, varargin)
% 
%  Properties:
%    - endo_list: List of endogenous variables.
%    - exo_list: List of exogenous variables.
%    - start_date: Last date of history.
%    - end_date: Last date of simulation.
%    - InitialRegime: Initial regime for forecasting (default is 1).
%    - Conditions: Conditioning information cell array {var_names, dates, values}.
% 
%  Dependent Properties:
%    - simul_periods: Number of periods for simulation.
%    - variable_list: List of exogenous and endogenous variables + regime.
% 
%  Methods:
%    - `append`: Add conditioning information to a pre-existing simulation plan.
%    - `export`: Compile the conditioning information in a `simplan` object for simulation or conditional forecasting.
% 
%  Description:
%    The `simplan` class is designed to facilitate the specification of
%    simulation and forecasting plans for Dynamic Stochastic General
%    Equilibrium (DSGE) models. It allows users to set conditions on
%    specific variables at certain dates, enabling detailed control over
%    the simulation process.
% 
%  Example:
%    % Create a `simplan` object
%    plan = simplan(model, horizon, initialRegime);
% 
%    % Add conditioning entries
%    plan = append(plan, 'inflation', someDate, 2.5);
%    plan = append(plan, 'interest_rate', someOtherDate, 1.8);
% 
%    % Export the compiled conditioning information
%    [histdb, np] = export(plan);
% 
%  See also: `append`, `export`
%
%    Documentation for simplan
%       doc simplan
%
%