%--- help for statespace/forecast ---
%
% FORECAST Forecast states and observations of state-space models
% 
%  Syntax:
% 
%    [Y,YMSE,X,XMSE] = forecast(Mdl,numPeriods,Y0)
%    [Y,YMSE,X,XMSE] = forecast(Mdl,numPeriods,name,value,...)
% 
%  Description:
% 
%    Generate multiple-period forecasts for observation vector y(t) and state 
%    vector x(t) for a general state-space model (SSM) of the form:
% 
%    State equation:       x(t) = A(t) * x(t-1) + B(t) * u(t)
%    Observation equation: y(t) = C(t) * x(t)   + D(t) * e(t)
% 
%    where u(t) and e(t) are uncorrelated, unit-variance white noise vector
%    processes. The length of x(t), y(t), u(t), and e(t) is m, n, k, and h, 
%    respectively.
% 
%    Additionally, forecast uncertainties are also generated.
% 
%  Input Arguments:
% 
%    Mdl - A state-space model, as created by the SSM constructor or 
%      SSM/ESTIMATE method.
% 
%    numPeriods - Positive, scalar, integer specifying the forecast horizon.
% 
%    Y0 - Observed response data to be forecasted. For time-invariant models 
%      in which the length of each observation vector (n) is the same, Y0 is 
%      a T-by-n matrix. For time-varying models in which the length of the 
%      observation vector changes, Y0 is a T-by-1 cell array in which each 
%      element contains a time-varying n-element vector of observations, y(t), 
%      associated with the corresponding period. The last observation is the 
%      most recent.
% 
%  Optional Input Name/Value Pairs:
% 
%    'A' Cell vector of forecasted state transition matrices in which each 
%        element is a matrix corresponding to a period in the forecast horizon 
%        at time t = 1,2,...,numPeriods. If the length of the state vector 
%        x(t) is constant, then each element is a square m-by-m matrix; 
%        however, if the length of x(t) changes, then some elements are 
%        non-square matrices. The length of the cell array must be at least 
%        numPeriods, and any elements beyond the forecast horizon are ignored.
%        By default, the last coefficient of the input SSM model Mdl is used 
%        for all future periods.
% 
%    'B' Cell vector of forecasted state disturbance loading matrices in which 
%        each element is a matrix corresponding to a period in the forecast 
%        horizon at time t = 1,2,...,numPeriods. If the lengths of the state 
%        vector x(t) and disturbance vector u(t) are constant, then each 
%        element is an m-by-k matrix; however, if the length of x(t) or u(t) 
%        changes, then the elements are matrices of various sizes. The length 
%        of the cell array must be at least numPeriods, and any elements 
%        beyond the forecast horizon are ignored. By default, the last 
%        coefficient of the input SSM model Mdl is used for all future periods.
% 
%    'C' Cell vector of forecasted measurement sensitivity matrices in which
%        each element is a matrix corresponding to a period in the forecast 
%        horizon at time t = 1,2,...,numPeriods. If the lengths of the 
%        observation vector y(t) and state vector x(t) are constant, then each 
%        element is an n-by-m matrix; however, if the length of y(t) or x(t) 
%        changes, then the elements are matrices of various sizes. The length 
%        of the cell array must be at least numPeriods, and any elements 
%        beyond the forecast horizon are ignored. By default, the last 
%        coefficient of the input SSM model Mdl is used for all future periods.
% 
%    'D' Cell vector of forecasted observation innovation matrices in which 
%        each element is a matrix corresponding to a period in the forecast 
%        horizon at time t = 1,2,...,numPeriods. If the lengths of the 
%        observation vector y(t) and innovation vector e(t) are constant, then 
%        each element is an n-by-h matrix; however, if the length of y(t) or 
%        e(t) changes, then the elements are matrices of various sizes. The 
%        length of the cell array must be at least numPeriods, and any 
%        elements beyond the forecast horizon are ignored. By default, the 
%        last coefficient of the input SSM model Mdl is used for all future 
%        periods.
% 
%    'Predictors0' T-by-d matrix of common predictor variables used to
%                  include a regression component in the observation equation. 
%                  Observations at time t are deflated such that
% 
%                  [y(t) - z(t)*b] = C * x(t) + D * e(t)
% 
%                  where z(t) is a vector of predictor variables and b is 
%                  the regression coefficient vector (see below). The default
%                  is an empty matrix (no regression component)
% 
%    'PredictorsF' numPeriods-by-d matrix of common predictor variables used to
%                  include a regression component in the observation equation. 
% 
%    'Beta'        d-by-n matrix of regression coefficients associated with
%                  predictors (see above). 
% 
%  Output Arguments:
% 
%    Y - Point forecasts of observations, E[y(t)|y(t-1),...,y(1)], for
%      t = 1,2,...,numPeriods. For time-invariant models in which the length 
%      of each observation vector (n) is the same, this is a numPeriods-by-n 
%      matrix. For time-varying models in which the length of the observation 
%      vector changes, this is a numPeriods-by-1 cell array in which each 
%      element contains a time-varying n-element vector of forecasts associated
%      with the corresponding period. 
% 
%    YMSE - Forecast error variances of future observations, 
%      Cov[y(t)|y(t-1),...,y(1)], for t = 1,2,...,numPeriods. For 
%      time-invariant models in which the length of each observation vector 
%      (n) is the same, this is a numPeriods-by-n matrix. For time-varying 
%      models in which the length of the observation vector changes, this is 
%      a numPeriods-by-1 cell array in which each element contains a 
%      time-varying n-element vector of forecast error variances associated 
%      with the corresponding period.
% 
%    X - Point forecasts of states, E[x(t)|y(t-1),...,y(1)], for t = 1,2,..., 
%      numPeriods. For time-invariant models in which the length of each state 
%      vector (m) is the same, this is a numPeriods-by-m matrix. For 
%      time-varying models in which the length of the state vector changes, 
%      this is a numPeriods-by-1 cell array in which each element contains a 
%      time-varying m-element vector of forecasts associated with the 
%      corresponding period. 
% 
%    XMSE - Forecast error variances of future states, 
%      Cov[x(t)|y(t-1),...,y(1)], for t = 1,2,...,numPeriods. For 
%      time-invariant models in which the length of each state vector (m) is 
%      the same, this is a numPeriods-by-m matrix. For time-varying models 
%      in which the length of the state vector changes, this is a 
%      numPeriods-by-1 cell array in which each element contains a time-varying 
%      m-element vector of forecast error variances associated with the 
%      corresponding period.
% 
%  See also SSM, FILTER, SMOOTH, SIMULATE, ESTIMATE.
%
%    Other functions named forecast
%
%       abstvar/forecast         garch/forecast
%       arima/forecast           generic/forecast
%       conjugateblm/forecast    gjr/forecast
%       customblm/forecast       regARIMA/forecast
%       diffuseblm/forecast      semiconjugateblm/forecast
%       egarch/forecast          varm/forecast
%       empiricalblm/forecast    vecm/forecast
%