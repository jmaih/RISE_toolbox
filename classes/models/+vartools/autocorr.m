% AUTOCORR Sample autocorrelation
% 
%  Syntax:
% 
%    [acf,lags] = autocorr(y)
%    ACFTbl = autocorr(Tbl)
%    [...,bounds] = autocorr(...)
%    [...,bounds,h] = autocorr(...)
%    [...] = autocorr(...,param,val,...)
%    [...] = autocorr(ax,...)
%    autocorr(...)
% 
%  Description:
% 
%    Compute the sample autocorrelation function (ACF) of a univariate 
%    time series. AUTOCORR optionally plots the ACF with confidence bounds.
% 
%  Input Arguments:
% 
%    y - Univariate time series data, specified as a numeric vector.
% 
%    Tbl - Time series data, specified as a table or timetable. Specify a
%        single series for y using the 'DataVariable' parameter.
% 
%    ax - Axes object in which to plot. If unspecified, AUTOCORR plots to
%        the current axes (gca).
% 
%    The function treats NaN values as missing completely at random.
% 
%  Optional Input Parameter Name/Value Arguments:
% 
%   'NumLags' Positive integer that determines the number of lags at which the
%             ACF is computed. The lags used to compute the ACF are 0:NumLags.
%             The default is min[20,N-1], where N is the effective sample size
%             of y.
% 
%   'NumMA'   For computing confidence bounds, a nonnegative integer less
%             than NumLags specifying the number of lags in a theoretical
%             MA(NumMA) model of y. For lags > NumMA, AUTOCORR uses
%             Bartlett's approximation [1] to compute the standard error
%             under the model assumption. The default is 0, in which case
%             the standard error is 1/sqrt(N), for Gaussian white noise.
% 
%   'NumSTD'  For computing confidence bounds, a nonnegative scalar multiple
%             specifying an interval of +/-(NumSTD) times the computed
%             standard error. The default is 2 (approximate 95% confidence).
% 
%   'DataVariable' Variable in Tbl to use for y, specified as a name in
%             Tbl.Properties.VariableNames. Variable names are character
%             vectors, string scalars, integers or logical vectors. The
%             default is the last variable in Tbl.
% 
%  Output Arguments:
% 
%    acf - Sample ACF. Vector of length NumLags+1 of values computed at lags
%        0,1,2,...,NumLags. For all y, acf(1) = 1 at lag 0.
% 
%    lags - Vector of lag numbers of length NumLags+1 used to compute acf.
% 
%    ACFTbl - When input is Tbl, outputs lags and acf are returned in
%        table ACFTbl.
% 
%    bounds - Two-element vector of approximate upper and lower confidence
%        bounds, assuming that y is an MA(NumMA) process.
% 
%    h - Vector of handles to plotted graphics objects. AUTOCORR plots the
%        ACF when the number of output arguments is 0 or 4.
% 
%  Notes:
% 
%   o If y is fully observed, without NaNs, AUTOCORR uses a Fourier
%     transform to compute the ACF in the frequency domain, then converts
%     back to the time domain using an inverse Fourier transform.
% 
%   o In the presence of NaNs, AUTOCORR computes the ACF at lag k in the
%     time domain, including in the sample average only those terms for
%     which the cross product y(t)*y(t+k) exists, so that the effective
%     sample size at any lag is a random variable.
% 
%  Example:
% 
%    % Create an MA(2) process from a sequence of 1000 Gaussian deviates,
%    % and assess whether the ACF is effectively zero for lags > 2:
% 
%    x = randn(1000,1);         % 1000 Gaussian deviates ~ N(0,1)
%    y = filter([1 -1 1],1,x);  % Create an MA(2) process
%    autocorr(y,'NumMA',2)      % Inspect the ACF with 95% confidence
% 
%  References:
% 
%    [1] Box, G. E. P., G. M. Jenkins, and G. C. Reinsel. Time Series
%        Analysis: Forecasting and Control. 3rd edition. Upper Saddle River,
%        NJ: Prentice-Hall, 1994.
% 
%    [2] Hamilton, J.D. Time Series Analysis. Princeton, NJ: Princeton
%        University Press, 1994.
% 
%  See also PARCORR, CROSSCORR, FILTER.
%
%    Documentation for autocorr
%       doc autocorr
%
%