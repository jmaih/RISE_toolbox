%--- help for abstvar/estimate ---
%
% --- help for abstvar/estimate ---
% 
%   Estimates VAR parameters using frequentist/Bayesian interpretation
%  
%   ::
%  
%      var = estimate(var);
%      var = estimate(var, data);
%      var = estimate(var, data, data_range);
%      var = estimate(var, data, data_range, prior);
%      var = estimate(var, data, data_range, prior, restrictions);
%      var = estimate(var, data, data_range, prior, restrictions,optimizer);
%      var = estimate(var, data, data_range, prior, restrictions,optimizer,is_fixed_regime);
%  
%   Args:
%  
%      var (var object): var object
%  
%      data (struct or ts obect): data for estimation
%  
%      data_range (serial): (optional) date_range
%  
%      prior (string): priors for parameters
%  
%         - If no prior is given, maximum likehood estimators (frequentist)
%         are set for parameters 
%         - With priors, posterior mode (Bayesian) values are set for parameters
%  
%      restrictions : restrictions. Refer to XXXXXXXX for analysis
%  
%      optimizer (char|function_handle|cell|{fmincon}) : optimization
%         procedure. Used with optimization is required. e.g. markov
%         switching. This can be the name of a standard matlab optimizer or
%         RISE optimization routine or a user-defined optimization procedure
%         available of the matlab search path. If the optimzer is provided as
%         a cell, then the first element of the cell is the name of the
%         optimizer or its handle and the remaining entries in the cell are
%         additional input arguments to the user-defined optimization
%         routine. A user-defined optimization function should have the
%         following syntax :: 
%  
%              [xfinal,ffinal,exitflag,H]=optimizer(fh,x0,lb,ub,options,varargin);
%  
%           That is, it accepts as inputs:
%  
%               - **fh**: the function to optimize
%               - **x0**: a vector column of initial values of the parameters
%               - **lb**: a vector column of lower bounds
%               - **ub**: a vector column of upper bounds
%               - **options**: a structure of options whose fields will be similar
%                 to matlab's optimset
%               - **varargin**: additional arguments to the user-defined
%                 optimization procedure
%  
%           That is, it provides as outputs:
%  
%               - **xfinal**: the vector of final values
%               - **ffinal**: the value of **fh** at **xfinal**
%               - **exitflag**: a flag similar to the ones provided by matlab's
%                 optimization functions.
%               - **H**: an estimate of the Hessian
%  
%      is_fixed_regime : (true|{false}): if true, the regimes are known in
%         advance. In that case the data should contain time series for a
%         variable called "hist_regimes"
%  
%   Returns:
%      : var object with parameters estimated based on data
%  
%   See also:
%      - :func:`posterior_mode <var.posterior_mode>`
%      - :func:`identification <var.identification>`
%  
% 
%     Other functions named estimate
% 
%        arima/estimate           gjr/estimate
%        conjugateblm/estimate    lassoblm/estimate
%        customblm/estimate       mixconjugateblm/estimate
%        diffuseblm/estimate      mixsemiconjugateblm/estimate
%        dsge/estimate            regARIMA/estimate
%        dssm/estimate            semiconjugateblm/estimate
%        egarch/estimate          ssm/estimate
%        empiricalblm/estimate    statespace/estimate
%        garch/estimate           varm/estimate
%        generic/estimate         vecm/estimate
% 
%    Other functions named estimate
%
%       arima/estimate           gjr/estimate
%       conjugateblm/estimate    lassoblm/estimate
%       customblm/estimate       mixconjugateblm/estimate
%       diffuseblm/estimate      mixsemiconjugateblm/estimate
%       dsge/estimate            regARIMA/estimate
%       dssm/estimate            semiconjugateblm/estimate
%       egarch/estimate          ssm/estimate
%       empiricalblm/estimate    statespace/estimate
%       garch/estimate           varm/estimate
%       generic/estimate         vecm/estimate
%