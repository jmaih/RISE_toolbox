%--- help for arima/simulate ---
%
% SIMULATE Simulate ARIMA model responses and conditional variances
% 
%  Syntax:
% 
%    [Y,E,V] = simulate(Mdl,numObs)
%    [Y,E,V] = simulate(Mdl,numObs,param1,val1,param2,val2,...)
% 
%  Description:
% 
%    Simulate sample paths of responses, innovations, and conditional 
%    variances of a univariate ARIMA process.
% 
%  Input Arguments:
% 
%    Mdl - ARIMA model specification object, as produced by the ARIMA 
%      constructor or ARIMA/ESTIMATE method.
% 
%    numObs - Positive integer indicating the number of observations (rows)
%      generated for each path of the outputs Y, E, and V.
% 
%  Optional Input Parameter Name/Value Pairs:
% 
%    'NumPaths'   Positive integer indicating the number of sample paths 
%                 (columns) generated for all simulated outputs. The default 
%                 is 1.
% 
%    'Y0'         Presample response data, providing initial values for the 
%                 model. Y0 is a column vector or a matrix. If Y0 is a 
%                 column vector, then it is applied to each simulated path. 
%                 If Y0 is a matrix, then it must have at least NumPaths 
%                 columns. Y0 may have any number of rows, provided at least 
%                 Mdl.P observations exist to initialize the model. If the
%                 number of rows exceeds Mdl.P, then only the most recent 
%                 Mdl.P observations are used. If the number of columns 
%                 exceeds NumPaths, then only the first NumPaths columns are 
%                 used. If Y0 is unspecified, any necessary presample 
%                 observations are set to the unconditional mean for stationary 
%                 AR processes, and to zero if the process is non-stationary 
%                 or contains a regression component. The last row contains 
%                 the most recent observation.
% 
%    'E0'         Mean-zero presample innovations, providing initial values 
%                 for the model. E0 is a column vector or a matrix. If E0 is 
%                 a column vector, then it is applied to each simulated path. 
%                 If E0 is a matrix, then it must have at least NumPaths
%                 columns. E0 may have any number of rows, provided sufficient 
%                 observations exist to initialize the ARIMA model as well 
%                 as any conditional variance model (the number of observations
%                 required is at least Mdl.Q, but may be more if a conditional 
%                 variance model is included). If the number of rows 
%                 exceeds the number necessary, then only the most recent 
%                 observations are used. If the number of columns exceeds 
%                 NumPaths, then only the first NumPaths columns are used. If
%                 no presample data is specified, any necessary observations 
%                 are set to zero. The last row contains the most recent 
%                 observation.
% 
%    'V0'         Positive presample conditional variances, providing initial
%                 values for any conditional variance model; if the variance 
%                 of the model is constant, then V0 is unnecessary. V0 is a 
%                 column vector or a matrix. If V0 is a column vector, then 
%                 it is applied to each simulated path. If V0 is a matrix, 
%                 then it must have at least NumPaths columns. V0 may have 
%                 any number of rows, provided sufficient observations exist 
%                 to initialize any conditional variance model. If the number 
%                 of rows exceeds the number necessary, then only the most 
%                 recent observations are used. If the number of columns 
%                 exceeds NumPaths, then only the first NumPaths columns are 
%                 used. If no presample variance data is specified, any 
%                 necessary observations are set to the unconditional variance 
%                 of the conditional variance process. The last row contains 
%                 the most recent observation.
% 
%    'X'          Matrix of predictor data used to include a regression 
%                 component in the conditional mean. Each column of X is a
%                 separate time series, and the last row of each contains
%                 the most recent observation of each series. The number of
%                 observations in X must equal or exceed numObs. When the 
%                 number of observations in X exceeds numObs, only the most 
%                 recent observations are used. If missing, the conditional 
%                 mean will have no regression component regardless of the 
%                 presence of any regression coefficients found in the model.
% 
%  Output Arguments:
% 
%    Y - numObs-by-NumPaths matrix of simulated response data.
% 
%    E - numObs-by-NumPaths matrix of simulated mean-zero innovations. 
% 
%    V - numObs-by-NumPaths matrix of conditional variances of the 
%      innovations in E.
% 
%  Notes:
% 
%    o Missing data values, indicated by NaNs, are removed from X by listwise 
%      deletion (i.e., any row in X with at least one NaN is removed), reducing
%      the effective sample size. The presample data Y0, E0, and V0 are merged
%      into a composite series, and any row of the combined series with at 
%      least one NaN is also removed by listwise deletion. The presample data 
%      is also synchronized such that the last (most recent) observation of 
%      each series occurs at the same time.
% 
%   o  Regression models included in the conditional mean are based on the
%      presence of the predictor matrix X. Although each column of the output 
%      time series represents a different path of the corresponding univariate 
%      stochastic process, the regression matrix X represents as a single 
%      path of a (possibly) multivariate time series matrix in which each 
%      column is a different time series. When the conditional mean has a 
%      regression component, the entire predictor matrix X is applied to 
%      every column of the output time series. 
% 
%  References:
% 
%    [1] Box, G. E. P., G. M. Jenkins, and G. C. Reinsel. Time Series
%        Analysis: Forecasting and Control. 3rd edition. Upper Saddle River,
%        NJ: Prentice-Hall, 1994.
% 
%    [2] Enders, W. Applied Econometric Time Series. Hoboken, NJ: John Wiley
%        & Sons, 1995.
% 
%    [3] Hamilton, J. D. Time Series Analysis. Princeton, NJ: Princeton
%        University Press, 1994.
% 
%  See also ARIMA, FORECAST, ESTIMATE, INFER, FILTER.
%
%    Reference page in Doc Center
%       doc arima/simulate
%
%    Other functions named simulate
%
%       conjugateblm/simulate    generic/simulate
%       customblm/simulate       gjr/simulate
%       diffuseblm/simulate      regARIMA/simulate
%       dsge/simulate            sde/simulate
%       dtmc/simulate            semiconjugateblm/simulate
%       egarch/simulate          ssm/simulate
%       empiricalblm/simulate    varm/simulate
%       garch/simulate           vecm/simulate
%