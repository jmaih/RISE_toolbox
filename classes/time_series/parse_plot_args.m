% PARSE_PLOT_ARGS Defines further plotting arguments for Rise time series objects.
% 
%  [this, rise_items, matlab_args] = parse_plot_args(varargin) is an
%  internal function that allows you to control various properties of the
%  plot for Rise time series objects.
% 
%  Input:
%    varargin (variable number of arguments) - Input arguments that define
%      the properties of the plot. These arguments can include:
%      - The range of values to be plotted: 'xrange', specified as a vector
%        or a date range in serial format. Default is an empty vector ([]).
%      - The number of tick marks on the plot: 'nticks', specified as a
%        positive integer. Default is determined by the plot_specs function.
%      - The vertical lines to be plotted: 'vline', specified as a string
%        representing a date or a range of dates. Default is an empty string ('').
%      - The horizontal lines to be plotted: 'hline', specified as a scalar
%        or a vector of numeric values. Default is an empty string ('').
%      - The log scale for the plot: 'logy', specified as a logical value.
%        Default is false.
%      - The date format for labeling the x-axis: 'date_format', specified as
%        a string. Default is an empty string ('').
% 
%  Output:
%    this - A time series object.
%    rise_items - A struct containing various properties for controlling the plot,
%      including 'xrange', 'nticks', 'vline', 'hline', 'logy', and 'date_format'.
%    matlab_args - Additional arguments to be passed to MATLAB's plotting functions.
% 
%  Notes:
%  - This function is intended for internal use and is called by other
%    functions in the Rise toolbox.
%  - The function allows you to customize the appearance of plots for Rise
%    time series objects, including specifying the range of values to be
%    plotted, the number of tick marks, the presence of vertical and
%    horizontal lines, and the use of a logarithmic scale.
%  - The 'xrange' argument can be used to specify the range of values to be
%    plotted on the x-axis. It can be specified as a vector of values or as a
%    date range in serial format. If only one value is provided, it is
%    interpreted as the end value and the start value is determined
%    automatically based on the frequency of the time series. If no 'xrange'
%    argument is provided, the entire range of the time series is plotted.
%  - The 'nticks' argument specifies the number of tick marks on the plot's
%    x-axis. It should be a positive integer. If not provided, the default
%    value is determined by the plot_specs function.
%  - The 'vline' argument is used to specify vertical lines to be plotted on
%    the plot. It can be specified as a string representing a single date or
%    a range of dates separated by commas. The dates should be in the same
%    frequency as the time series being plotted.
%  - The 'hline' argument is used to specify horizontal lines to be plotted
%    on the plot. It can be specified as a scalar or a vector of numeric
%    values.
%  - The 'logy' argument specifies whether to use a logarithmic scale for
%    the y-axis. It should be a logical value (true or false). If true, the
%    logarithm of the time series values will be taken before plotting. Note
%    that if any data values are less than or equal to zero, a warning will
%    be displayed.
%  - The 'date_format' argument is used to specify the format of the dates
%    displayed on the x-axis. It should be a string in the format expected by
%    MATLAB's datestr function. If not provided, the default format is used.
% 
%  Example:
%    data = randn(100, 1);
%    dates = rq(1990,1);
%    db = ts(dates, data);
%    [this, rise_items, matlab_args] = parse_plot_args(db, 'nticks', 10, 'logy', true);
% 
%  See also: PLOT_SPECS, RISE_TIME_SERIES, DATE2SERIAL, DATESTR
%