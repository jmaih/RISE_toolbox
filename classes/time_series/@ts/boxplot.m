function h=boxplot(varargin)
% H1 line
%
% ::
%
%
% Args:
%
% Returns:
%    :
%
% Note:
%
% Example:
%
%    See also:


%  boxplot Displays box plots of multiple time series in a frame.
%     boxplot(db) produces a box plot of the data in db.  If db has many
%     variables, there is one box per column, and if db is univariate there
%     is just one box. On each box, the central mark is the median, the
%     edges of the box are the 25th and 75th percentiles, the whiskers
%     extend to the most extreme datapoints the algorithm considers to be
%     not outliers, and the outliers are plotted individually.  
%     
%     boxplot(xrange,db) selects the data matching the time represented by
%     xrange
%
%     boxplot(..., 'PARAM1', val1, 'PARAM2', val2, ...) specifies optional
%     parameter name/value pairs.
%       'plotstyle'     'traditional' (default), or 'compact' to specify a
%                       box style designed for plots with many groups.  The
%                       plotstyle changes the defaults for some other
%                       parameters, as described below.
%  
%       'boxstyle'      'outline' (default) to draw an unfilled box with
%                       dashed lines for whiskers, or 'filled' to draw a
%                       narrow filled box with solid lines for whiskers.
%       'colorgroup'    One or more grouping variables, of the same type as 
%                       permitted for G, specifying that the box color should
%                       change when the specified variables change.  Default
%                       is [] for no box color change.
%       'colors'        Colors for boxes, specified as a single color (such
%                       as 'r' or [1 0 0]) or multiple colors (such as 'rgbm'
%                       or a three-column matrix of RGB values).  The sequence
%                       is replicated or truncated as required, so for example
%                       'rb' gives boxes that alternate in color.  Default
%                       when no 'colorgroup' is specified is to use the same
%                       color scheme for all boxes.  Default with
%                       'colorgroup' is a modified hsv colormap. 
%       'datalim'       A two-element vector containing lower and upper limits,
%                       used by 'extrememode' to determine which points are
%                       extreme.  Default is [-Inf Inf].
%       'extrememode'   'clip' (default) to move data outside the 'datalim'
%                       limits to the limit, or 'compress' to distribute such
%                       points evenly in a region just outside the limit,
%                       retaining the relative order of the points.  A
%                       dotted line marks the limit if any points are outside
%                       it, and two gray lines mark the compression region if
%                       any points are compressed.  Values at +/-Inf can be
%                       clipped or compressed, but NaNs still do not appear
%                       on the plot.  Box notches are drawn to scale and may
%                       extend beyond the bounds if the median is inside the
%                       limit; they are not drawn if the median is outside
%                       the limits.  
%       'factordirection' 'data' (default) to arrange the factors with the
%                       first value next to the origin, 'list' to arrange the
%                       factors left-to-right if on the x axis or top-to-
%                       bottom if on the y axis, or 'auto' to use 'data' for
%                       numeric grouping variables and 'list' for strings.
%       'fullfactors'   'off' (default) to have one group for each unique row
%                       of G, or 'on' to create a group for each possible 
%                       combination of group variable values, including
%                       combinations that do not appear in the data.
%       'factorseparator' Specifies which factors should have their values 
%                       separated by a grid line.  The value may be 'auto' or
%                       a vector of grouping variable numbers.  For example,
%                       [1 2] adds a separator line when the first or second
%                       grouping variable changes value.  'auto' is [] for
%                       one grouping variable and [1] for two or more
%                       grouping variables. Default is [].
%       'factorgap'     Specifies an extra gap to leave between boxes when
%                       the corresponding grouping factor changes value,
%                       expressed as a percentage of the width of the plot.
%                       For example, with [3 1], the gap is 3% of the width
%                       of the plot between groups with different values of
%                       the first grouping variable, and 1% between groups
%                       with the same value of the first grouping variable
%                       but different values for the second.  'auto'
%                       specifies that boxplot should choose a gap
%                       automatically.  Default is [].
%       'grouporder'    Order of groups for plotting, specified as a cell
%                       array of strings.  With multiple grouping variables,
%                       separate values within each string with a comma.
%                       Using categorical arrays as grouping variables is an
%                       easier way to control the order of the boxes.
%       'jitter'        Maximum distance D to displace outliers along the
%                       factor axis by a uniform random amount, in order to
%                       make duplicate points visible.  D = 1 makes the
%                       jitter regions just touch between the closest
%                       adjacent groups.  The default is 0.
%       'labels'        Character array, cell array of strings, or numeric
%                       vector of box labels.  May have one label per group
%                       or per db value.  Multiple label variables may be
%                       specified via a numeric matrix or a cell array
%                       containing any of these types.
%       'labelorientation' 'horizontal' (default) for horizontal labels, or
%                       'inline' to draw the labels vertically when
%                       'orientation' has its default 'vertical' value.
%       'labelverbosity'  'all' (default) to display every label, 'minor' to
%                       display a label for a factor only when that factor
%                       has a different value from the previous group, or
%                       'majorminor' to display a label for a factor when
%                       that factor or any factor major to it has a
%                       different value from the previous group.
%       'medianstyle'   'line' (default) to draw a line for the median, or
%                       'target' to draw a black dot inside a white circle.
%       'notch'         'on' to draw comparison intervals using notches
%                       ('plotstyle' is 'traditional) or triangular markers
%                       ('plotstyle' is 'compact'), 'marker' to draw them
%                       using triangular markers, or 'off' (default) to omit
%                       them.  Two medians are significantly different at the
%                       5% level if their intervals do not overlap.  The
%                       interval endpoints are the extremes of the notches or
%                       the centers of the triangular markers.  When the
%                       sample size is small, notches may extend beyond the
%                       end of the box.
%       'orientation'   'vertical' (default) to plot db on the y axis, or
%                       'horizontal' to plot db on the x axis.
%       'outliersize'   Size of marker used for outliers, in points.
%                       Default is 6.
%       'positions'     Box positions specified as a numeric vector with one
%                       entry per group or db value (default 1:NGROUPS when the
%                       number of groups is NGROUPS).
%       'symbol'        Symbol and color to use for outliers, using the same 
%                       values as the LineSpec parameter S in PLOT.  Default
%                       is 'r+'. If the symbol is omitted then the outliers
%                       are invisible; if the color is omitted then the
%                       outliers have the same color as their corresponding
%                       box.  Any line specification in S is ignored.
%       'whisker'       Maximum whisker length W.  Default is W=1.5.  Points
%                       are drawn as outliers if they are larger than
%                       Q3+W*(Q3-Q1) or smaller than Q1-W*(Q3-Q1), where Q1
%                       and Q3 are the 25th and 75th percentiles, respectively.
%                       The default value 1.5 corresponds to approximately +/-
%                       2.7 sigma and 99.3 coverage if the data are normally
%                       distributed.  The plotted whisker extends to the
%                       adjacent value, which is the most extreme data value
%                       that is not an outlier. Set 'whisker' to 0 to give no
%                       whiskers and to make every point outside of Q1 and Q3
%                       an outlier.
%       'widths'        A scalar or vector of box widths to use when the
%                       'boxstyle' is 'outline'.  The default is half of the
%                       minimum separation between boxes, which is .5 when
%                       the 'positions' argument takes its default value.
%                       The list of values is replicated or truncated as
%                       necessary.
%  
%     When the 'plotstyle' parameter takes the value 'compact', then the
%     default values for other parameters are the following:
%         boxstyle - 'filled'            labelverbosity - 'majorminor'
%         factorgap - 'auto'             medianstyle - 'target'
%         factorseparator - 'auto'       outliersize - 4
%         jitter - 0.5                   symbol - 'o'
%         labelorientation - 'inline'        
%  
%     You can see the data values and group names by using the data cursor
%     tool, available from the figure window.  The data cursor shows the
%     original values of any points affected by the 'datalim' parameter.  You
%     can label the specific group to which an outlier belongs using the gname
%     function.
%  
%     To modify the properties of box components, use findobj using tags to
%     find their handles as in one of the examples below.  The tag names
%     depend on the plotstyle and are:
%  
%        all styles:  'Box', 'Outliers'
%        traditional: 'Median', 'Upper Whisker', 'Lower Whisker',
%                     'Upper Adjacent Value', 'Lower Adjacent Value', 
%        compact:     'Whisker', 'MedianOuter', 'MedianInner'
%        when 'notch' is 'marker':
%                     'NotchLo', 'NotchHi'
%     In addition to those matlab properties, RISE adds further properties,
%     which allow to control for. See parse_plot_args
                   
tmp=utils.plot.myplot(@boxplot,varargin{:});

if nargout
    h=tmp;
end

end
