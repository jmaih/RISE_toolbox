% HPFILTER Hodrick-Prescott filter for trend and cyclical components
% 
%  Syntax:
% 
% 	[Trend,Cyclical] = hpfilter(Y)
% 	[Trend,Cyclical] = hpfilter(Y,smoothing)
% 	hpfilter(...)
% 
%  Description:
% 
%    Separate one or more time series into trend and cyclical components
%    with the Hodrick-Prescott filter. If no output arguments are specified,
%    HPFILTER displays a plot of the series and trend (with cycles removed).
%    The plot can be used to help select a smoothing parameter.
% 
%  Input Arguments:
% 
% 	Y - Time series data. Y may be a vector or a matrix. If Y is a vector,
% 	  it represents a single series. If Y is a numObs-by-numSeries matrix,
% 	  it represents numObs observations of numSeries series, with
% 	  observations across any row assumed to occur at the same time. The
% 	  last observation of any series is assumed to be the most recent.
% 
%  Optional Input Argument:
% 
% 	smoothing - Either a scalar to be applied to all series or a vector of
% 	    length numSeries with values to be applied to corresponding series.
% 	    The default is 1600, which is suggested	in [1] for quarterly data.
% 	    If smoothing is 0, no smoothing occurs.	As the smoothing parameter
% 	    increases, the smoothed series approaches a straight line. If
% 	    smoothing is Inf, the series is detrended.
% 
%  Output Arguments:
% 
% 	Trend - Trend component of Y, the same size as Y.
% 	Cyclical - Cyclical component of Y, the same size as Y.
% 
%  Notes:
% 
% 	o The Hodrick-Prescott filter separates a time series into a trend
% 	  component and a cyclical component such that
% 
% 		Y = Trend + Cyclical
% 
% 	  The filter is equivalent to a cubic spline smoother, where the
% 	  smoothed portion is in Trend.
% 
% 	o [1] suggests values for the smoothing parameter that depend upon
% 	  the periodicity of the data:
% 
% 		Periodicity     smoothing
%        -----------     ---------
% 		Yearly			100
% 		Quarterly		1600
% 		Monthly			14400
% 
%    o The Hodrick-Prescott filter can produce anomalous endpoint effects in
%      very high-frequency data and should never be used for extrapolation.
% 
%  Reference:
% 
% 	[1] Hodrick, R. J., and E. C. Prescott. "Postwar U.S. Business Cycles:
% 		An Empirical Investigation." Journal of Money, Credit, and Banking.
% 		Vol. 29, 1997, pp. 1-16.
%
%    Reference page in Doc Center
%       doc hpfilter
%
%    Other functions named hpfilter
%
%       ts/hpfilter
%