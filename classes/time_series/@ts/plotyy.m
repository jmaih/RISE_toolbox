function varargout=plotyy(varargin)
% Plot figure with y tick labels on the left and right
%
% ::
%
%    [ax, h1, h2] = plotyy(y1,y2);
%    [ax, h1, h2] = plotyy(xrange, y1, y2);
%    [ax, h1, h2] = plotyy(xrange, y1, y2, func);
%    [ax, h1, h2] = plotyy(xrange, y1, y2, func1, func2);
%
% Args:
%    y1 (ts object): time series to plot with left axis tick
%    y2 (ts object): time series to plot with right axis tick
%    xrange: date range to plot. See ts for formats
%    func (string/func-handle): use the function to make the plots instead
%    func1 (string/func-handle): use the function to make the plot for y1
%    func2 (string/func-handle): use the function to make the plot for y2
%
% Returns:
%    :
%
%       - ax: axis handle
%       - h1: handle to the graphics object corresponding to y1
%       - h2: handle to the graphics object corresponding to y2
%
% Example:
%    ::
%
%       [ax,h1,h2]=plotyy('1994m7:1997m1',...
%           db(:,'v1',1),...
%           db(:,'v2',1),...
%           'nticks',10,...
%           'date_format',17)
%           xrotate(90)
%
% Note:
%     In addition to those matlab properties, RISE adds further properties,
%     which allow to control for. See parse_plot_args
%


[varargout{1:nargout}]=utils.plot.myplot(@plotyy,varargin{:});

end