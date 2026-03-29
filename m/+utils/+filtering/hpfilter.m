% HPFILTER Hodrick-Prescott filter for trend and cyclical components
% 
%  Syntax:
% 
%    [Trend,Cyclical] = hpfilter(Y)
%    [TTbl,CTbl] = hpfilter(Tbl)
%    [...,h] = hpfilter(...)
%    [...] = hpfilter(...,param,val,...)
%    [...] = hpfilter(ax,...)
%    hpfilter(...)
% 
%  Description:
% 
%    Separate time series into additive trend and cyclical components with
%    the Hodrick-Prescott filter [1]. HPFILTER optionally plots the series
%    together with the smoothed (cycles removed) trend components.
% 
%  Input Arguments:
% 
%    Y - Time series data, specified as a numObs-by-numVars numeric matrix.
% 
%    Tbl - Time series data, specified as a table or timetable. Specify
%        series for Y using the 'DataVariables' parameter.
% 
%    ax - Axes object in which to plot. If unspecified, HPFILTER plots to
%        the current axes (gca).
% 
%    The function removes NaN values indicating missing observations.
% 
%  Optional Input Parameter Name/Value Arguments:
% 
%    NAME            VALUE
% 
%    'Smoothing'     Smoothing parameter, specified as a nonnegative scalar
%                    or vector. Scalar values are applied to all series in
%                    Y. Vector values of length numVars are applied to
%                    corresponding series in Y. As smoothing increases, the
%                    trend approaches a straight line. The default is 1600,
%                    suggested in both [1] and [3] for quarterly data.
% 
%    'FilterType'    Type of filter, specified as a string or character
%                    vector. Values are 'two-sided' [1] or 'one-sided' [4].
%                    The default is 'two-sided'.
% 
%    'DataVariables' Variables in Tbl to use for Y, specified as names in
%                    Tbl.Properties.VariableNames. Variable names are cell
%                    vectors of character vectors, string vectors, integer
%                    vectors or logical vectors. The default is all
%                    variables in Tbl.
% 
%  Output Arguments:
% 
%    Trend - Smoothed trend component of Y, the same size as Y.
% 
%    Cyclical - Cyclical component of Y, the same size as Y.
% 
%    TTbl - When input is Tbl, output Trend is returned in tabular TTbl, the
%        same type as Tbl.
% 
%    CTbl - When input is Tbl, output Cyclical is returned in tabular CTbl,
%        the same type as Tbl.
% 
%    h - Vector of handles to plotted graphics objects. HPFILTER plots the
%        data and trend when the number of output arguments is 0 or 3.
% 
%  Notes:
% 
%    o The Hodrick-Prescott filter separates time series Y into a trend
%      component and a cyclical component such that
% 
%        Y = Trend + Cyclical
% 
%      The method implements a high-pass filter for the cycle that penalizes
%      variations in the trend to a degree determined by the smoothing
%      parameter [1].
% 
%    o Hodrick and Prescott [1] suggest values for the smoothing parameter
%      that depend upon the periodicity of the data:
% 
%        Periodicity     Smoothing
%        -----------     ---------
%        Yearly			100
%        Quarterly		1600
%        Monthly			14400
% 
%      Ravn and Uhlig [3] suggest adjustments to these values:
% 
%        Periodicity     Smoothing
%        -----------     ---------
%        Yearly			6.25
%        Quarterly		1600
%        Monthly			129600
% 
%      In practice, a vector of smoothing parameters allows for testing of
%      alternatives. The plot produced by HPFILTER is useful for comparison
%      of results.
% 
%    o The default two-sided filter uses future values of the input series
%      to compute outputs at time t and is typically applied to historical
%      data. It may produce anomalous end effects unsuitable for forecasting
%      [4]. The one-sided filter, by contrast, is causal, using only current
%      and previous values of the input series. As a result, the one-sided
%      filter does not revise outputs when new data becomes available.
% 
%  Example:
% 
%    % Filter annual GNP data:
% 
%    load Data_NelsonPlosser
%    [TTbl,CTbl,h] = hpfilter(DataTable,'DataVariables',["GNPR" "GNPN"],...
%                             'Smoothing',6.25);
% 
%  References:
% 
%    [1] Hodrick, R. J., and E. C. Prescott. "Postwar U.S. Business Cycles:
%        An Empirical Investigation." Journal of Money, Credit, and Banking.
%        Vol. 29, 1997, pp. 1-16.
% 
%    [2] Hodrick, R.J. "An Exploration of Trend-Cycle Decomposition
%        Methodologies in Simulated Data." National Bureau of Economic
%        Research. Working Paper 26750. 2020.
% 
%    [3] Ravn, M. O. and H. Uhlig. "On Adjusting the Hodrick-Prescott Filter
%        for the Frequency of Observations." The Review of Economics and
%        Statistics. Vol. 84, No. 2, 2002, pp. 371-376.
% 
%    [4] Stock, J. and M. Watson. "Forecasting inflation." Journal of
%        Monetary Economics. Vol. 44, No. 2, 1999, pp. 293-335.
%  
%  See also BKFILTER, CFFILTER, HFILTER.
%
%    Documentation for hpfilter
%       doc hpfilter
%
%    Other uses of hpfilter
%
%       ts/hpfilter
%