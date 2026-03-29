%--- help for dssm/fevd ---
%
% FEVD Forecast error variance decomposition of state-space models
% 
%  Syntax:
% 
%    For fully specified models:
%    Decomposition = fevd(Mdl)
%    Decomposition = fevd(Mdl,name,value,...)
% 
%    For partially specified models:
%    [...] = fevd(Mdl,'Params',estParams,'EstParamCov',EstParamCov,...)
%    [...] = fevd(Mdl,'Params',estParams,'EstParamCov',EstParamCov,name,value,...)
%    where estParams and EstParamCov are returned by ESTIMATE method.
% 
%  Description:
% 
%    For observation vector y(t) and state vector x(t), a general 
%    state-space model (SSM) takes the form:
% 
%    State equation:       x(t) = A(t) * x(t-1) + B(t) * u(t)
%    Observation equation: y(t) = C(t) * x(t)   + D(t) * e(t)
% 
%    where u(t) and e(t) are uncorrelated, unit-variance white noise vector
%    processes. The lengths of x(t), y(t), u(t), and e(t) are m, n, k, and h, 
%    respectively.
% 
%    The forecast error variance decomposition (FEVD) attributes the
%    volatility of y(t) to component-wise shocks to u(t). FEVD provides
%    information about the relative importance of each shock in affecting
%    the forecast error variance of y(t).
% 
%  Input Arguments:
% 
%    Mdl        - A state-space model (@ssm or @dssm object).
% 
%  Optional Input Name/Value Pairs:
% 
%    'Params'      Estimates of the unknown parameters when Mdl is partially
%                  specified, specified a numeric vector. For an explicitly
%                  created model, arrange the elements of Params to
%                  correspond to hits of a column-wise search of NaNs in the
%                  state-space model parameters in this order: A, B, C, D,
%                  Mean0, and Cov0. For an implicitly created model,
%                  ParamMap is solely responsible for mapping the parameter
%                  vector into model parameters. This input is only
%                  necessary if the model has unknown parameters.
% 
%    'EstParamCov' Estimated covariance matrix of unknown parameters when 
%                  Mdl is partially specified, specified as a positive
%                  semidefinite numeric matrix. ESTIMATE returns the
%                  estimated parameter covariance matrix of Mdl in the
%                  appropriate form. The default is empty.
% 
%    'NumPeriods'  Number of periods in the FEVD, specified as a positive 
%                  integer. The default is 20.
% 
%    'NumPaths'    Number of sample paths (trials) to generate to estimate 
%                  confidence bounds, specified as a positive integer. 
%                  The default is 1000.
% 
%    'Confidence'  Confidence level, specified as a nonnegative scalar. 
%                  Confidence (C) is such that 100*(1-C)/2 percent of the
%                  variance decompositions lie below and above the lower and
%                  upper confidence bounds, respectively. Confidence must be
%                  between 0 and 1. The default is 0.95.
% 
%  Output Arguments:
% 
%    Decomposition - numPeriods-by-k-by-n array of variance decomposition.
%                 The FEVD in element (t,i,j) is the contribution to the
%                 variance decomposition of variable j attributable to a
%                 shock to variable i at time t for t = 1,2,...,numPeriods.
% 
%    Lower -      numPeriods-by-k-by-n array of pointwise lower bounds of
%                 the decomposition.
% 
%    Upper -      numPeriods-by-k-by-n array of pointwise upper bounds of
%                 the decomposition.
%                  
% 
%  Notes:
% 
%  o Forecast variance of y(t) is caused by noise in the transition
%    and measurement equations. For a non-zero coefficient matrix D, the
%    decomposition does not sum up one (along the second dimension). The
%    remaining portion is attributable to D*D'.
% 
%  o Models with time-varying coefficients have a time-varying FEVD.
%    The timing convention is that the state shock occurs in period 1.
%    FEVD is calculated for period 1, 2, and so on.
% 
%  See also FILTER, SMOOTH, ESTIMATE, FORECAST.
%
%    Documentation for dssm/fevd
%       doc dssm/fevd
%
%    Other uses of fevd
%
%       ssm/fevd    statespace/fevd    varm/fevd    vecm/fevd
%