%--- help for statespace/refine ---
%
% REFINE Refine initial parameters to aid estimation of state-space models
% 
%  Syntax:
% 
%    Output = refine(Mdl,Y,params0)
%    Output = refine(Mdl,Y,params0,name,value,...)
%    refine(...)
% 
%  Description:
% 
%    For observation vector y(t), state vector x(t), and uncorrelated, unit-
%    variance white noise vector processes u(t) and e(t), maximum likelihood 
%    parameter estimation of the following state-space model:
% 
%    State equation:       x(t) = A(t) * x(t-1) + B(t) * u(t)
%    Observation equation: y(t) = C(t) * x(t)   + D(t) * e(t)
% 
%    is often sensitive to the initial parameters specified by the user.
% 
%    In the event a preliminary estimation fails to converge, or converges to
%    an unsatisfactory solution, this utility function refines the original 
%    initial parameters in an attempt to improve model estimation results.
% 
%    In the model above, the length of x(t), y(t), u(t), and e(t) is m, n, k, 
%    and h, respectively.
% 
%  Input Arguments:
% 
%    Mdl - A state-space model with unknown parameters to estimate, created 
%      by the SSM constructor.
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
%  Output Argument:
% 
%    Output - Output structure array in which each element has the following 
%      fields:
% 
%      o Description    A brief description of the refinement algorithm.
%                       Current algorithms are:
% 
%                       'Quasi-Newton'
%                       'Nelder-Mead simplex'
%                       'Loose bound interior point'
%                       'Starting value perturbation'
%                       'Starting value shrinkage'
% 
%      o Parameters:    Vector of initial parameter values.
% 
%      o LogLikelihood: Log-likelihood associated with the initial vector.
% 
%  Notes:
% 
%   o When called with no output argument, a tabular display of summary
%     information related to the resulting vector of initial parameter values
%     associated with each algorithm is printed to the screen. If the Output
%     structure array is requested, then no summary display is printed.
% 
%   o The various refined initial parameter vectors may look similar to each 
%     other and to the original (see params0 input above). Manually select a
%     candidate vector of initial parameter values that makes economic sense
%     and is associated with a relatively large log-likelihood value.
% 
%   o If a refinement attempt is unsuccessful, its corresponding log-likelihood 
%     is set to -Inf and its parameter vector is empty. Error messages are 
%     displayed on the screen.
%  
%  See also SSM, DSSM, ESTIMATE, FILTER, SMOOTH, FORECAST, SIMULATE, SIMSMOOTH.
%