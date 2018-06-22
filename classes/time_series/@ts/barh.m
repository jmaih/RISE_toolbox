function plot_handle=barh(varargin)
% Make bar graph of the time series
%
% ::
%
%    varargout = bar(X,db,varargin);
%
% Args:
%    - X : Date description
%
%       - string of dates, e.g., X='1990:2000'
%       - serial dates, e.g., X=date2serial('1990Q1'):date2serial('2000Q1')
%
%    - db (ts object): Time series object with data
%    - options: Options need to come in pairs
%
%       - 'figsize' ([value, value]): the figure size (multiple plots) with default [3,3]
%       - 'figtitle' (string): the figure title (multiple plots) with default ''
%       - the number of tick marks (integer): 'nticks', with default 8
%       - 'date_format': the date format (see matlab's datestr) with default ''
%       - 'logy' (bool): the log scale with default false
%       - 'secondary_y': the list of variables going to the secondary y axis
%       - 'subplots' (bool): the flag for multiple plots with default false
%
% Notes:
%
%   - bar(db) uses the dates within db.  The colors are set by the colormap.
%   - bar(X,db,WIDTH) or bar(db,WIDTH) specifies the width of the bars. Values
%     of WIDTH > 1, produce overlapped bars.  The default value is WIDTH=0.8
%   - bar(...,'grouped') produces the default vertical grouped bar chart.
%     bar(...,'stacked') produces a vertical stacked bar chart.
%     bar(...,LINESPEC) uses the line color specified (one of 'rgbymckw').
%   - H = bar(...) returns a vector of handles to barseries objects.
%   - Use SHADING FACETED to put edges on the bars.  Use SHADING FLAT to
%     turn them off.
%
% Example:
%    ::
%
%       bar('1994m7:1997m1',db(:,:,1),...
%           'figsize',[2,2],...
%           'figtitle','no title',...
%           'nticks',10,...
%           'legend',{'v1','v2'},...
%           'legend_loc','BO',...
%           'logy',true,...
%           'secondary_y',{'v1','v4'},...
%           'subplots',true,...
%           'linewidth',2,...
%           'date_format',17);
%       xrotate(90)
%
% See also:
%    - hist
%    - plot
%    - barh
%    - bar3
%    - bar3h
%

plot_handle=utils.plot.myplot(@barh,varargin{:});

end