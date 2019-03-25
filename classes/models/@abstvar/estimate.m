%--- help for statespace/estimate ---
%
% ESTIMATE Maximum likelihood parameter estimation of state-space models 
% 
%  Syntax:
% 
%    [EstMdl,estParams,EstParamCov,logL,Output] = estimate(Mdl,Y,params0)
%    [EstMdl,estParams,EstParamCov,logL,Output] = estimate(Mdl,Y,params0,
%                                                          name,value,...)
% 
%  Description:
% 
%    For observation vector y(t) and state vector x(t), estimate parameters 
%    of the following general state-space model by maximum likelihood:
% 
%    State equation:       x(t) = A(t) * x(t-1) + B(t) * u(t)
%    Observation equation: y(t) = C(t) * x(t)   + D(t) * e(t)
% 
%    where u(t) and e(t) are uncorrelated, unit-variance white noise vector
%    processes. The length of x(t), y(t), u(t), and e(t) is m, n, k, and h, 
%    respectively.
% 
%    For models created explicitly by specifying coefficients A, B, C, and D,
%    unknown parameters to estimate are identified by the presence of NaNs in
%    these model coefficients. Unknown parameters identified by the presence 
%    of NaNs in the mean vector (Mean0) and covariance matrix (Cov0) of initial
%    states x(0) are optionally estimated as well.
% 
%    For models created implicitly by specifying a parameter mapping function 
%    ParamMap, the mapping function is responsible for managing the presence 
%    and placement of unknown parameters. In the implicit approach, the mapping
%    function alone defines the model, and is particularly convenient for 
%    estimating complex models and for imposing certain parameter constraints. 
%    Moreover, in more general settings in which the initial states are also 
%    determined by unknown parameters, ParamMap may include additional output 
%    arguments; refer to the SSM/DSSM constructor for more details.
% 
%  Input Arguments:
% 
%    Mdl - A state-space model with unknown parameters to estimate, created 
%      by the SSM/DSSM constructor.
% 
%    Y - Observed response data to which the model is fit. For time-invariant 
%      models in which the length of each observation vector (n) is the same, 
%      Y is a T-by-n matrix. For time-varying models in which the length of 
%      the observation vector changes, Y is a T-by-1 cell array in which 
%      each element contains a time-varying n-element vector of observations, 
%      y(t), associated with the corresponding period. The last observation 
%      is the most recent.
% 
%    params0 - A vector containing the initial values of unknown 
%      parameters associated with model coefficients A, B, C, and D, and 
%      optionally the mean vector (Mean0) and covariance matrix (Cov0) of 
%      initial states x(0), estimated by maximum likelihood. For models created 
%      explicitly, parameters mapped to NaN values are found by a column-wise 
%      search of A, followed by B, then C, then D, and finally Mean0 and Cov0. 
%      For models created implicitly, the parameter function ParamMap is 
%      solely responsible for mapping the initial parameter vector into model 
%      coefficients A, B, C, and D, as well as additional information 
%      regarding initial states and types if necessary. 
% 
%  Optional Input Name/Value Pairs:
% 
%    'Univariate'  Logical value indicating whether to use the univariate 
%                  treatment of a multivariate series. The default is false.
% 
%    'SquareRoot'  Logical value indicating whether to use the square-root 
%                  filter. The default is false.
%                  SSM only. No effects on DSSM.
% 
%    'Tolerance'   A small, non-negative variance tolerance that controls 
%                  whether an observed series is ignored if its forecast
%                  uncertainty falls below this threshold. Setting tolerance 
%                  to a small number, say 1e-15, may help overcome numerical 
%                  problems of the filter. The default is 0.
% 
%    'Predictors'  T-by-d matrix of common predictor variables used to
%                  include a regression component in the observation equation. 
%                  Observations at time t are deflated such that
% 
%                  [y(t) - z(t)*b] = C * x(t) + D * e(t)
% 
%                  where z(t) is a vector of predictor variables and b is 
%                  the regression coefficient vector (see below). The default
%                  is an empty matrix (no regression component)
% 
%    'Beta0'       d-by-n matrix of initial values of regression coefficients 
%                  associated with predictors (see above). If the model contains
%                  a regression component, coefficients are estimated along 
%                  with other unknown parameters in A, B, C, D, Mean0, and 
%                  Cov0; the default initial values are obtained by ordinary 
%                  least squares (OLS) by regressing Y on the explanatory 
%                  variables.
% 
%    'SwitchTime'  a positive integer that specifies the date after which the
%                  diffuse filter switches to the standard filter.
%                  The default is the earliest date when the smoothed
%                  initial states have a full-rank covariance matrix.
%                  DSSM only. No effects on SSM.
% 
%    'CovMethod'   String or character vector indicating the method for
%                  computing the asymptotic parameter error covariance
%                  matrix of estimated parameters. Values are:
% 
%                  VALUE             METHOD
% 
%                  o 'hessian'       Negative inverted Hessian matrix
%                  o 'opg'           Outer product of gradients (default)
%                  o 'sandwich'      Both Hessian and outer product of gradients
% 
%    'Options'     Optimization options created with OPTIMOPTIONS. If 
%                  specified, default optimization parameters are replaced 
%                  by those in options. The default OPTIMOPTIONS object depends
%                  on the optimization function used to estimate parameters.
%                  For constrained optimization, FMINCON is called and the 
%                  algorithm used is interior point; for unconstrained 
%                  optimization, FMINUNC is called and the algorithm is 
%                  quasi-Newton. See documentation for OPTIMOPTIONS, FMINCON, 
%                  and FMINUNC for details.
% 
%    'Display'    String vector or cell vector of character vectors
%                 indicating what information to display in the command
%                 window. Values are:
%   
%                 VALUE           DISPLAY
%   
%                 o 'off'         No display to the command window. 
%   
%                 o 'params'      Display maximum likelihood parameter 
%                                 estimates, standard errors, and t statistics.
%                                 This is the default.
% 
%                 o 'iter'        Display iterative optimization information.
% 
%                 o 'diagnostics' Display optimization diagnostics.
% 
%                 o 'full'        Display 'params', 'iter', and 'diagnostics'.            
% 
%    The following optional name/value pairs are related to constrained
%    optimization performed by FMINCON (see FMINCON for details):
% 
%    'Aineq'  Linear inequality matrix, such that for solution vector X 
%             Aineq*X <= bineq. The number of rows is determined by the 
%             number of constraints, and the number of columns by the number 
%             of estimated parameters. Columns are ordered by estimated 
%             parameters in A, B, C, D, Mean0, Cov0, and finally regression 
%             coefficients for models with a regression component.
% 
%    'bineq'  Linear inequality column vector, such that for solution vector
%             X Aineq*X <= bineq. 
% 
%    'Aeq'    Linear equality matrix, such that for solution vector X 
%             Aeq*X = beq. The number of rows is determined by the 
%             number of constraints, and the number of columns by the number 
%             of estimated parameters. Columns are ordered by estimated 
%             parameters in A, B, C, D, Mean0, Cov0, and finally regression 
%             coefficients for models with a regression component.
% 
%    'beq'    Linear equality column vector, such that for solution vector X 
%             Aeq*X = beq. 
% 
%    'lb'     Lower bounds column vector on estimated parameters. Vector 
%             elements are ordered by estimated parameters in A, B, C, D, 
%             Mean0, Cov0, and finally regression coefficients for models with 
%             a regression component.
% 
%    'ub'     Upper bounds column vector on estimated parameters. Vector 
%             elements are ordered by estimated parameters in A, B, C, D, 
%             Mean0, Cov0, and finally regression coefficients for models with 
%             a regression component.
% 
%  Output Arguments:
% 
%    EstMdl - A fitted state-space model with estimated model coefficients. 
%      Provided the optimization converged successfully, the model is explicit 
%      in its coefficients A, B, C, D, Mean0, and Cov0 regardless of whether 
%      the input model (Mdl) was created explicitly or implicitly.
% 
%    estParams - Vector of estimated parameter. The elements are ordered such
%      that any estimated parameters in A appear first, then B, then C, then
%      D, and finally Mean0 and Cov0. Additionally, if the observation equation
%      includes a regression component, then estimates of any regression 
%      coefficients appear last.
% 
%    EstParamCov - Variance-covariance matrix of estimated parameters. The
%      rows/columns are ordered such that variances/covariances of any 
%      estimated parameters in A appear first, then B, then C, then D, and 
%      finally Mean0 and Cov0. Additionally, if the observation equation
%      includes a regression component, then variances/covariances of any 
%      estimated regression coefficients appear last.
% 
%    logL - Log-likelihood of the observations.
% 
%    Output - Output structure with the following fields:
% 
%       o ExitFlag - Optimization exit flag that describes the exit condition.
%           See FMINCON or FMINUNC for additional details.
% 
%       o Options - Optimization options, created with OPTIMOPTIONS, and used
%           by the optimizer.
% 
%  Notes:
% 
%    o Missing observations in Y are indicated by NaNs.
% 
%    o Although creating a model explicitly by directly specifying parameters
%      (A, B, C, D, etc.) with NaN placeholders to indicate parameters to 
%      estimate is more convenient than specifying a user-defined mapping 
%      function ParamMap, the utility of an explicit approach is limited in 
%      that each estimated parameter affects and is uniquely associated with 
%      a single element of a coefficient matrix. 
% 
%    o The option of univariate treatment requires a diagonal D(t)*D(t)'.
% 
%    o The state-space model itself does not store regression data or 
%      coefficients. Instead, a regression component is used to deflate the 
%      observations, such that the deflated data is given by y(t) - z(t)*b. 
% 
%    o Regression coefficients in the observation equation are estimated along 
%      with other unknown parameters, and so parameter constraints may be 
%      placed on regression coefficients as well other parameters by specifying 
%      appropriate entries in 'Aineq', 'bineq', etc.
% 
%    o If observations are multivariate, then the same regressors apply to all.
% 
%    o If the state equation requires predictors, or each individual component 
%      of the observation equation requires a set of distinct predictors, try
%      one of the following methods:
% 
%        o Expand the states by a constant one.
%        o Expand the states by predictors.
%        o Return 8 output arguments from the user-supplied function 
%          ParamMap, the last of which is the deflated observation data.
% 
%    o Regression components in the observation equation are allowed only for
%      time-invariant observations, in which the input observation series Y
%      is of constant length and specified as a T-by-n matrix. Y specified as
%      a T-by-1 cell array indicates a time-varying observation series whose 
%      length may change over time. If the length of each observation y(t)
%      changes, then it is unclear which regression coefficients are needed
%      to deflate a given observation.
% 
%    o If no constraints are specified, FMINUNC is used; otherwise, FMINCON 
%      is used for constrained optimization. Therefore, when specifying the 
%      input Options, the input should be consistent with the solver (see 
%      OPTIMOPTIONS, FMINCON, and FMINUNC for details).
% 
%      Whenever possible, it is recommended to avoid equality/inequality 
%      constraints by reparameterizing the model. For example, various 
%      parameters may be exponentiated using EXP to ensure positivity.
% 
%    o The observations in the starting periods are treated as if they were
%      presample data for exact initialization of the diffuse initial
%      states. We define the likelihood function as the joint density
%      of the observations after the algorithm is switched to
%      the standard Kalman filter. The variable 'SwitchTime' determines when
%      the standard Kalman filter starts.
%   
%  See also SSM, DSSM, FILTER, SMOOTH, FORECAST, SIMULATE, SIMSMOOTH.
%
%    Other functions named estimate
%
%       abstvar/estimate         garch/estimate
%       arima/estimate           generic/estimate
%       conjugateblm/estimate    gjr/estimate
%       customblm/estimate       regARIMA/estimate
%       diffuseblm/estimate      semiconjugateblm/estimate
%       dsge/estimate            ssm/estimate
%       dssm/estimate            varm/estimate
%       egarch/estimate          vecm/estimate
%       empiricalblm/estimate
%