%--- help for statespace/irf ---
%
% IRF Impulse response function of state-space models
% 
%  Syntax:
% 
%    For fully specified models:
%    [ResponseY,ResponseX] = irf(Mdl)
%    [...] = irf(Mdl,name,value,...)
% 
%    For partially specified models:
%    [...] = irf(Mdl,'Params',estParams,'EstParamCov',EstParamCov,...)
%    [...] = irf(Mdl,'Params',estParams,'EstParamCov',EstParamCov,name,value,...)
%    where estParams and EstParamCov are returned by ESTIMATE.
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
%    processes. The length of x(t), y(t), u(t), and e(t) is m, n, k, and h, 
%    respectively.
% 
%    Impulse response functions (IRF) measure the dynamic effects on the
%    state and measurement series x(t), y(t), t = 1,2,..., when an
%    unanticipated change to u(1) occurs.
%    
% 
%  Input Arguments:
% 
%    Mdl        - A state-space model, as created by the model constructor
%                 or ESTIMATE.
% 
%  Optional Input Name/Value Pairs:
% 
%    'Params'      A vector of parameter values associated with unknown
%                  parameters found in model coefficients A, B, C, and D, and
%                  optionally the mean vector (Mean0) and covariance matrix 
%                  (Cov0) of initial states x(0), estimated by maximum 
%                  likelihood. For models created explicitly, parameters mapped
%                  to NaN values are found by a column-wise search of A, 
%                  followed by B, then C, then D, and finally Mean0 and Cov0. 
%                  For models created implicitly, the parameter mapping function 
%                  ParamMap is solely responsible for mapping the initial 
%                  parameter vector into model coefficients A, B, C, and D, 
%                  as well as additional information regarding initial states 
%                  and types if necessary. This input is only necessary if the
%                  model has unknown parameters. 
% 
%    'EstParamCov' The output EstParamCov (estimator covariance matrix)
%                 returned by ESTIMATE for computing lower and upper
%                 bounds. The default is empty.
% 
%    'NumPeriods' Positive integer indicating the number of periods in the
%                 IRF. The default is 20.
% 
%    'NumPaths'   Number of sample paths (trials) to generate for computing 
%                 lower and upper bounds. The default is 1000.
% 
%    'Confidence' Positive scalar for the confidence level of the bounds. 
%                 The default is 0.95.
% 
%    'Cumulative' Logical value indicating cumulative responses instead of
%                 period-by-period responses. The default is false.
% 
%    'Method'     String that specifies the irf estimation algorithm.
%                 o 'repeated-multiplication' (default).
%                 o 'eigendecomposition': the software attempts to use the
%                   eigenvalues of the time-invariant transition matrix for
%                   computing impulse response functions.
% 
%  Output Argument:
% 
%    ResponseY -  numPeriods-by-k-by-n array of dynamic responses of Y.
%                 If a unit of shock to variable i occurs in period 1,
%                 then the response of variable j in period t is given by
%                 the element (t,i,j).
% 
%    ResponseX -  numPeriods-by-k-by-m array of dynamic responses of X.
%                 If a unit of shock to variable i occurs in period 1,
%                 then the response of variable j in period t is given by
%                 the element (t,i,j).
% 
%    LowerY -     numPeriods-by-k-by-n array of pointwise lower bounds 
%                 of dynamic responses of Y. 
% 
%    UpperY -     numPeriods-by-k-by-n array of pointwise upper bounds 
%                 of dynamic responses of Y.
% 
%    LowerX -     numPeriods-by-k-by-m array of pointwise lower bounds 
%                 of dynamic responses of X. 
% 
%    UpperX -     numPeriods-by-k-by-m array of pointwise upper bounds 
%                 of dynamic responses of X. 
% 
%  Notes:
% 
%  o The timing convention is that the state shock occurs in period 1.
%    The responses are calculated for period 1, 2, and so on.
% 
%  o Models with time-varying coefficients have a time-varying IRF.
% 
%  o If 'EstParamCov' (estimator covariance matrix) is not available,
%    lower and upper bounds overlap.
% 
%  See also FILTER, SMOOTH, ESTIMATE, FORECAST.
%
%    Other uses of irf
%
%       dsge/irf    generic/irf    varm/irf        vecm/irf
%       dssm/irf    ssm/irf        varsimul/irf
%