%--- help for dsge/counterfactual ---
%
%  counterfactual Computes counterfactual history of a nonlinear DSGE model
% 
%  Syntax
%  -------
%  ::
% 
%    counterf=counterfactual(m)
%    counterf=counterfactual(m,sim_engine)
%    counterf=counterfactual(m,sim_engine,nsim)
%    counterf=counterfactual(m,sim_engine,nsim,shock_names)
%    counterf=counterfactual(m,sim_engine,nsim,shock_names,varargin)
% 
%  Inputs
%  -------
% 
%  - m : [rise|dsge] model(s) for which to compute the
%    counterfactual. m could be a vector of models. The computation of the
%    counterfactual requires a computation of smoothed history. Hence either
%    m should contain all the information needed for that purpose or the
%    information should be passed along through varargin below
% 
%  - sim_engine : [empty|function handle] function computing forecasts
% 
%  - nsim : [empty|numeric|{100}] Number of simulation to consider for the
%    integration exercise. nsim is automatically set to 1 if the model is
%    detected to be solved at order 1 and not contain regime switches.
% 
%  - shock_names : [empty|char|cellstr] list of shocks to consider in the
%    computation of the counterfactual
% 
%  - varargin : additional information needed for the computation of the
%    smoothed quantities through filtration.
% 
%  Outputs
%  --------
% 
%  - counterf : [struct|cell array] structure or cell array of structures
%    with the counterfactuals in each model and each solution.
% 
%  Remarks
%  --------
% 
%  - **N.B** : For a nonlinear model (e.g linear but switching), the type of
%    decomposition/counterfactual we do for linear/linearized
%    constant-parameter models is not feasible. RISE uses a monte carlo
%    integration to provide an approximation to the
%    decomposition/counterfactual
% 
%  - **N.B** : if m is a vector of models, then each model should return a
%    unique solution (implying a unique filtration), else the concatenation
%    of counterfactuals will fail. In that case it is better to run one
%    model at a time.
% 
%  - **N.B** : When a model has multiple solutions, there is no guarantee
%    that every solution will have succeed at filtration. For this reason,
%    the result of each counterfactual is stored in a separate cell. Empty
%    cells represent the solutions for which the filtration could not be
%    obtained.
% 
%  Examples
%  ---------
% 
%  See also: dsge/historical_decomposition
%